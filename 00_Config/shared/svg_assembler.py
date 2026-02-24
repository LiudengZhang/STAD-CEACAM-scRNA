#!/usr/bin/env python3
"""
VectorAssembler — Compose SVG panels into vector PDF figures.

Uses direct SVG XML composition via lxml for the merged output, and
cairosvg (Python) for PDF and PNG conversion. This approach handles
SVGs of any size without the svglib parsing bottleneck or rsvg-convert
element limits.

Output:
  - SVG with all panels as nested vector elements (Illustrator-editable)
  - PDF converted from merged SVG via cairosvg (vector)
  - PNG preview at 300 DPI
"""

import base64
import os
import re
import tempfile
from io import BytesIO
from pathlib import Path

import cairocffi as cairo
import cairosvg
from lxml import etree
from PIL import Image


SVGNS = 'http://www.w3.org/2000/svg'
XLINKNS = 'http://www.w3.org/1999/xlink'
INKSCAPENS = 'http://www.inkscape.org/namespaces/inkscape'


class VectorAssembler:
    """Assemble SVG/PNG panels into a fully vector figure.

    Coordinate system: origin at top-left, units in mm.
    All panel placements use (x_mm, y_mm) from top-left corner.

    SVG panels are embedded by direct XML nesting — no svglib parsing
    required, so even large (50+ MB) SVGs are handled in seconds.
    PNG panels are embedded as base64 <image> elements.
    """

    RENDER_DPI = 300

    def __init__(self, width_mm, height_mm, title=None):
        self.width_mm = width_mm
        self.height_mm = height_mm
        self.panels_placed = []
        self._panel_sources = []
        self._elements = []  # List of (type, data) tuples
        self._panel_counter = 0
        self._tmp_dir = tempfile.mkdtemp(prefix='vec_asm_')

        if title:
            self.add_text(title, x_mm=width_mm / 2, y_mm=2,
                          font_size=10, anchor='middle', font_weight='bold')

    def place_panel(self, label, svg_path, x_mm, y_mm, w_mm, h_mm,
                    label_offset_x=1, label_offset_y=3, force_png=False,
                    force_vector=False):
        """Place an SVG (or PNG fallback) panel at the specified position.

        Args:
            label: Panel label (e.g. 'A'). Pass None to skip label.
            svg_path: Path to SVG file. If not found, tries .png sibling.
            x_mm, y_mm: Top-left corner position in mm.
            w_mm, h_mm: Target width and height in mm.
            label_offset_x: Label x offset from panel left edge (mm).
            label_offset_y: Label y offset from panel top edge (mm).
            force_png: If True, skip SVG and use PNG directly.
            force_vector: If True, render as vector in hybrid compositing
                          regardless of file size.
        """
        svg_path = Path(svg_path)

        if svg_path.suffix == '.png':
            png_path = svg_path
            svg_path = svg_path.with_suffix('.svg')
        else:
            png_path = svg_path.with_suffix('.png')

        if not force_png and svg_path.exists():
            self._panel_counter += 1
            self._elements.append(('svg', {
                'path': str(svg_path),
                'x_mm': x_mm, 'y_mm': y_mm,
                'w_mm': w_mm, 'h_mm': h_mm,
                'panel_id': self._panel_counter,
                'force_vector': force_vector,
                'label': label,
            }))
            source = 'SVG(vector)'
        elif png_path.exists():
            self._elements.append(('png', {
                'path': str(png_path),
                'x_mm': x_mm, 'y_mm': y_mm,
                'w_mm': w_mm, 'h_mm': h_mm,
                'label': label,
            }))
            source = 'PNG'
        else:
            print(f"  WARNING: Neither {svg_path.name} nor {png_path.name} found")
            return

        if label:
            lx = x_mm + label_offset_x
            ly = y_mm + label_offset_y
            self.add_text(label, lx, ly, font_size=8, font_weight='normal')

        self.panels_placed.append(label or svg_path.stem)
        self._panel_sources.append((label or svg_path.stem, source))
        print(f"  Panel {label or svg_path.stem}: {source} ({w_mm:.1f} x {h_mm:.1f} mm)")

    def add_text(self, text, x_mm, y_mm, font_size=8, font_name='Liberation Sans',
                 font_weight='normal', anchor='start', color=None, rotation=0):
        """Add vector text at the specified position.

        Args:
            text: Text string.
            x_mm, y_mm: Position in mm from top-left.
            font_size: Font size in points.
            font_name: Base font name.
            font_weight: 'normal' or 'bold'.
            anchor: 'start', 'middle', or 'end'.
            color: Color string or reportlab Color object.
            rotation: Rotation in degrees (counter-clockwise).
        """
        color_str = self._color_to_svg(color)
        self._elements.append(('text', {
            'text': text,
            'x_mm': x_mm, 'y_mm': y_mm,
            'font_size': font_size,
            'font_name': font_name,
            'font_weight': font_weight,
            'anchor': anchor,
            'color': color_str,
            'rotation': rotation,
        }))

    @staticmethod
    def _color_to_svg(color):
        """Convert a color value to SVG color string."""
        if color is None:
            return 'black'
        if isinstance(color, str):
            return color
        # Handle reportlab Color objects
        try:
            r = int(color.red * 255)
            g = int(color.green * 255)
            b = int(color.blue * 255)
            return f'#{r:02x}{g:02x}{b:02x}'
        except (AttributeError, TypeError):
            return 'black'

    def _get_svg_viewbox(self, panel_root):
        """Extract viewBox from SVG root, computing from width/height if needed."""
        vb = panel_root.get('viewBox', '')
        if vb:
            return vb

        w = panel_root.get('width', '')
        h = panel_root.get('height', '')
        if w and h:
            w_val = re.match(r'([\d.]+)', w.strip())
            h_val = re.match(r'([\d.]+)', h.strip())
            if w_val and h_val:
                return f'0 0 {w_val.group(1)} {h_val.group(1)}'
        return ''

    def _namespace_ids(self, elem, prefix):
        """Prefix all id= attributes and url(#...) / href='#...' references
        to prevent ID conflicts between panels in the merged SVG."""
        id_pattern = re.compile(r'url\(#([^)]+)\)')

        for node in elem.iter():
            # Rename id attributes
            node_id = node.get('id')
            if node_id is not None:
                node.set('id', f'{prefix}_{node_id}')

            # Update references in all attributes
            for attr_name in list(node.attrib):
                val = node.attrib[attr_name]
                if val is None:
                    continue
                # url(#id) references in clip-path, mask, fill, etc.
                if 'url(#' in val:
                    node.attrib[attr_name] = id_pattern.sub(
                        lambda m: f'url(#{prefix}_{m.group(1)})', val)
                # xlink:href="#id" and href="#id" references
                if attr_name in ('{http://www.w3.org/1999/xlink}href', 'href'):
                    if val.startswith('#'):
                        node.attrib[attr_name] = f'#{prefix}_{val[1:]}'

            # Update url(#...) in inline style attribute
            style = node.get('style') or ''
            if 'url(#' in style:
                node.set('style', id_pattern.sub(
                    lambda m: f'url(#{prefix}_{m.group(1)})', style))

    def _build_svg(self, output_path):
        """Build merged SVG by direct XML composition."""
        nsmap = {None: SVGNS, 'xlink': XLINKNS, 'inkscape': INKSCAPENS}
        root = etree.Element('{%s}svg' % SVGNS, nsmap=nsmap)
        root.set('width', f'{self.width_mm}mm')
        root.set('height', f'{self.height_mm}mm')
        root.set('viewBox', f'0 0 {self.width_mm} {self.height_mm}')
        root.set('version', '1.1')

        # White background
        bg = etree.SubElement(root, '{%s}rect' % SVGNS)
        bg.set('width', '100%')
        bg.set('height', '100%')
        bg.set('fill', 'white')

        for etype, data in self._elements:
            if etype == 'svg':
                self._embed_svg_panel(root, data)
            elif etype == 'png':
                self._embed_png_panel(root, data)
            elif etype == 'text':
                self._embed_text(root, data)

        tree = etree.ElementTree(root)
        tree.write(str(output_path), xml_declaration=True,
                   encoding='UTF-8', pretty_print=False)

    def _embed_svg_panel(self, root, data):
        """Embed an SVG panel as a nested <svg> element with ID namespacing."""
        svg_path = data['path']
        x, y, w, h = data['x_mm'], data['y_mm'], data['w_mm'], data['h_mm']
        panel_id = data['panel_id']
        label = data.get('label')
        prefix = f'p{panel_id}'

        parser = etree.XMLParser(recover=True, huge_tree=True)
        try:
            tree = etree.parse(svg_path, parser)
        except Exception as e:
            print(f"  WARNING: lxml failed to parse {Path(svg_path).name}: {e}")
            png_path = Path(svg_path).with_suffix('.png')
            if png_path.exists():
                self._embed_png_panel(root, {**data, 'path': str(png_path)})
            return

        panel_root = tree.getroot()
        vb = self._get_svg_viewbox(panel_root)

        # Wrap in <g> layer for Illustrator/Inkscape layer support
        layer_id = f'Panel_{label}' if label else f'Panel_{panel_id}'
        layer_label = f'Panel {label}' if label else f'Panel {panel_id}'
        group = etree.SubElement(root, '{%s}g' % SVGNS)
        group.set('id', layer_id)
        group.set('{%s}groupmode' % INKSCAPENS, 'layer')
        group.set('{%s}label' % INKSCAPENS, layer_label)

        # Create nested <svg> at target position and size
        nested = etree.SubElement(group, '{%s}svg' % SVGNS)
        nested.set('x', f'{x}')
        nested.set('y', f'{y}')
        nested.set('width', f'{w}')
        nested.set('height', f'{h}')
        if vb:
            nested.set('viewBox', vb)
        nested.set('preserveAspectRatio', 'xMidYMid meet')

        # Copy all children from the panel SVG
        for child in list(panel_root):
            nested.append(child)

        # Namespace all IDs to prevent conflicts between panels
        self._namespace_ids(nested, prefix)

    def _embed_png_panel(self, root, data):
        """Embed a PNG as a base64 <image> element."""
        png_path = data['path']
        x, y, w, h = data['x_mm'], data['y_mm'], data['w_mm'], data['h_mm']
        label = data.get('label')
        panel_id = data.get('panel_id', id(data))

        with open(png_path, 'rb') as f:
            b64 = base64.b64encode(f.read()).decode('ascii')

        # Wrap in <g> layer for Illustrator/Inkscape layer support
        layer_id = f'Panel_{label}' if label else f'Panel_{panel_id}'
        layer_label = f'Panel {label}' if label else f'Panel {panel_id}'
        group = etree.SubElement(root, '{%s}g' % SVGNS)
        group.set('id', layer_id)
        group.set('{%s}groupmode' % INKSCAPENS, 'layer')
        group.set('{%s}label' % INKSCAPENS, layer_label)

        img = etree.SubElement(group, '{%s}image' % SVGNS)
        img.set('x', f'{x}')
        img.set('y', f'{y}')
        img.set('width', f'{w}')
        img.set('height', f'{h}')
        img.set('preserveAspectRatio', 'xMidYMid meet')
        img.set('{%s}href' % XLINKNS, f'data:image/png;base64,{b64}')

    def _embed_text(self, root, data):
        """Embed text as an SVG <text> element."""
        x, y = data['x_mm'], data['y_mm']
        font_size = data['font_size']
        font_name = data['font_name']
        font_weight = data['font_weight']
        anchor = data.get('anchor', 'start')
        color = data.get('color', 'black')
        rotation = data.get('rotation', 0)

        anchor_map = {'start': 'start', 'middle': 'middle', 'end': 'end'}
        svg_anchor = anchor_map.get(anchor, 'start')

        # Convert pt to mm (1pt = 0.3528mm) for the mm-based viewBox
        font_size_mm = font_size * 0.3528

        txt = etree.SubElement(root, '{%s}text' % SVGNS)
        txt.set('x', f'{x}')
        txt.set('y', f'{y}')
        txt.set('font-family', f'Liberation Sans, Arial, Helvetica, sans-serif')
        txt.set('font-size', f'{font_size_mm}')
        if font_weight == 'bold':
            txt.set('font-weight', 'bold')
        txt.set('text-anchor', svg_anchor)
        txt.set('fill', color)

        if rotation:
            txt.set('transform', f'rotate({-rotation}, {x}, {y})')

        txt.text = data['text']

    # Merged SVGs below this threshold use single-pass cairosvg (fastest)
    _SMALL_SVG_THRESHOLD = 20_000_000  # 20 MB
    # Individual panel SVGs below this threshold render as vector in compositing
    _VECTOR_PANEL_THRESHOLD = 5_000_000  # 5 MB

    def save(self, stem, output_dir=None, svg_out=True):
        """Save the assembled figure as SVG, PDF, and PNG preview.

        Uses direct SVG composition (lxml) for the merged SVG output.
        For PDF/PNG conversion:
          - Small figures (<20MB): cairosvg on the merged SVG (vector PDF).
          - Large figures: hybrid panel-by-panel compositing — small panels
            render as vector (RecordingSurface), large panels as 300 DPI
            raster. Text/labels always vector.

        Args:
            stem: Output filename stem (e.g. 'Figure_04_merged').
            output_dir: Directory for output. Defaults to current directory.
            svg_out: If True, also output a combined SVG (for Illustrator).

        Returns:
            Tuple of (pdf_path, png_path).
        """
        if output_dir is None:
            output_dir = Path('.')
        else:
            output_dir = Path(output_dir)

        svg_path = output_dir / f'{stem}.svg'
        pdf_path = output_dir / f'{stem}.pdf'
        png_path = output_dir / f'{stem}.png'

        # Build merged SVG directly via lxml (handles any SVG size)
        self._build_svg(svg_path)
        svg_size = svg_path.stat().st_size
        print(f"\nSaved combined SVG: {svg_path} ({svg_size / 1e6:.1f} MB)")

        if svg_size < self._SMALL_SVG_THRESHOLD:
            # Small figure: cairosvg on merged SVG (all-vector PDF)
            self._save_via_cairosvg(svg_path, pdf_path, png_path)
        else:
            # Large figure: hybrid compositing (vector small panels + raster large panels)
            print(f"  Large SVG ({svg_size / 1e6:.0f} MB) — using hybrid compositing")
            self._save_hybrid(pdf_path, png_path)

        # Report
        print(f"\nCanvas: {self.width_mm:.1f} x {self.height_mm:.1f} mm")
        print(f"Panels placed: {len(self.panels_placed)}")

        n_vector = sum(1 for _, s in self._panel_sources if 'vector' in s)
        n_raster = len(self._panel_sources) - n_vector
        print(f"Vector panels: {n_vector} | Raster panels: {n_raster}")
        if n_raster > 0:
            raster_labels = [lbl for lbl, s in self._panel_sources if 'vector' not in s]
            print(f"  Raster panels: {', '.join(raster_labels)}")

        # Cleanup temp files
        import shutil
        if os.path.exists(self._tmp_dir):
            shutil.rmtree(self._tmp_dir, ignore_errors=True)

        return str(pdf_path), str(png_path)

    def _save_via_cairosvg(self, svg_path, pdf_path, png_path):
        """Convert merged SVG to PDF and PNG using cairosvg (vector output)."""
        try:
            cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
            print(f"Saved vector PDF: {pdf_path}")
        except Exception as e:
            print(f"  cairosvg PDF failed: {e}")

        px_per_mm = self.RENDER_DPI / 25.4
        target_w_px = int(self.width_mm * px_per_mm)
        try:
            cairosvg.svg2png(url=str(svg_path), write_to=str(png_path),
                             output_width=target_w_px, background_color='white')
            print(f"Saved PNG preview: {png_path}")
        except Exception as e:
            print(f"  cairosvg PNG failed: {e}")

    def _render_svg_to_recording(self, svg_path):
        """Render an SVG to a cairocffi RecordingSurface (vector).

        Uses cairosvg's internal Tree + Surface with a subclassed surface
        that outputs to a RecordingSurface instead of a file.

        Returns:
            (RecordingSurface, width_pt, height_pt)
        """
        from cairosvg.surface import PDFSurface as _CairoSVGPDF
        from cairosvg.surface import Tree as _Tree

        class _RecSurface(_CairoSVGPDF):
            def _create_surface(self, width, height):
                rec = cairo.RecordingSurface(cairo.CONTENT_COLOR_ALPHA, None)
                return rec, width, height

        tree = _Tree(url=str(svg_path))
        surface = _RecSurface(tree, None, dpi=96)
        return surface.cairo, surface.width, surface.height

    def _save_hybrid(self, pdf_path, png_path):
        """Hybrid panel-by-panel compositing: vector for small panels, raster for large.

        SVG panels smaller than _VECTOR_PANEL_THRESHOLD are rendered as vector
        via RecordingSurface. Larger SVGs and PNG panels are rasterized at
        300 DPI. Text/labels are always vector on the PDF.
        """
        mm_to_pt = 72.0 / 25.4
        w_pt = self.width_mm * mm_to_pt
        h_pt = self.height_mm * mm_to_pt

        px_per_mm = self.RENDER_DPI / 25.4
        canvas_w = int(self.width_mm * px_per_mm)
        canvas_h = int(self.height_mm * px_per_mm)

        # Create PDF surface via cairocffi
        pdf_surface = cairo.PDFSurface(str(pdf_path), w_pt, h_pt)
        pdf_ctx = cairo.Context(pdf_surface)

        # White background
        pdf_ctx.set_source_rgb(1, 1, 1)
        pdf_ctx.rectangle(0, 0, w_pt, h_pt)
        pdf_ctx.fill()

        # PNG canvas via Pillow
        png_canvas = Image.new('RGB', (canvas_w, canvas_h), 'white')

        for etype, data in self._elements:
            if etype in ('svg', 'png'):
                x_mm, y_mm = data['x_mm'], data['y_mm']
                w_mm, h_mm = data['w_mm'], data['h_mm']
                panel_path = data['path']

                x_pt = x_mm * mm_to_pt
                y_pt = y_mm * mm_to_pt
                w_target_pt = w_mm * mm_to_pt
                h_target_pt = h_mm * mm_to_pt

                w_px = max(1, int(w_mm * px_per_mm))
                h_px = max(1, int(h_mm * px_per_mm))

                # Decide: vector or raster?
                use_vector = (
                    etype == 'svg'
                    and (data.get('force_vector', False)
                         or os.path.getsize(panel_path) < self._VECTOR_PANEL_THRESHOLD)
                )

                if use_vector:
                    # --- Vector rendering via RecordingSurface ---
                    recording, rec_w, rec_h = self._render_svg_to_recording(
                        panel_path)
                    scale_x = w_target_pt / rec_w if rec_w else 1
                    scale_y = h_target_pt / rec_h if rec_h else 1

                    pdf_ctx.save()
                    pdf_ctx.translate(x_pt, y_pt)
                    pdf_ctx.scale(scale_x, scale_y)
                    pdf_ctx.set_source_surface(recording, 0, 0)
                    pdf_ctx.paint()
                    pdf_ctx.restore()

                    panel_name = Path(panel_path).stem
                    print(f"    {panel_name}: vector")

                    # PNG preview: render via cairosvg
                    png_bytes = cairosvg.svg2png(
                        url=panel_path, output_width=w_px,
                        background_color='white')
                else:
                    # --- Raster rendering at 300 DPI ---
                    if etype == 'svg':
                        png_bytes = cairosvg.svg2png(
                            url=panel_path, output_width=w_px,
                            background_color='white')
                    else:
                        with Image.open(panel_path) as img:
                            img_resized = img.convert('RGB').resize(
                                (w_px, h_px), Image.LANCZOS)
                            buf = BytesIO()
                            img_resized.save(buf, 'PNG')
                            png_bytes = buf.getvalue()

                    tmp_png = os.path.join(self._tmp_dir,
                                           f'panel_{id(data)}.png')
                    with open(tmp_png, 'wb') as f:
                        f.write(png_bytes)
                    img_surface = cairo.ImageSurface.create_from_png(tmp_png)

                    scale_x = w_target_pt / img_surface.get_width()
                    scale_y = h_target_pt / img_surface.get_height()

                    pdf_ctx.save()
                    pdf_ctx.translate(x_pt, y_pt)
                    pdf_ctx.scale(scale_x, scale_y)
                    pdf_ctx.set_source_surface(img_surface)
                    pdf_ctx.paint()
                    pdf_ctx.restore()

                    os.unlink(tmp_png)

                    panel_name = Path(panel_path).stem
                    print(f"    {panel_name}: raster (300 DPI)")

                # --- Draw onto PNG canvas ---
                panel_img = Image.open(BytesIO(png_bytes))
                panel_img = panel_img.convert('RGB').resize(
                    (w_px, h_px), Image.LANCZOS)
                x_px_pos = int(x_mm * px_per_mm)
                y_px_pos = int(y_mm * px_per_mm)
                png_canvas.paste(panel_img, (x_px_pos, y_px_pos))

            elif etype == 'text':
                # Draw text on PDF surface (always vector)
                x_pt = data['x_mm'] * mm_to_pt
                y_pt = data['y_mm'] * mm_to_pt
                font_size_pt = data['font_size']
                pdf_ctx.save()
                pdf_ctx.set_source_rgb(0, 0, 0)
                pdf_ctx.select_font_face(
                    'Liberation Sans',
                    cairo.FONT_SLANT_NORMAL,
                    cairo.FONT_WEIGHT_BOLD if data['font_weight'] == 'bold'
                    else cairo.FONT_WEIGHT_NORMAL)
                pdf_ctx.set_font_size(font_size_pt)
                pdf_ctx.move_to(x_pt, y_pt + font_size_pt)
                pdf_ctx.show_text(data['text'])
                pdf_ctx.restore()

                # Draw text on PNG canvas
                from PIL import ImageDraw, ImageFont
                draw = ImageDraw.Draw(png_canvas)
                x_px = int(data['x_mm'] * px_per_mm)
                y_px = int(data['y_mm'] * px_per_mm)
                font_size_px = int(data['font_size'] * 0.3528 * px_per_mm)
                try:
                    font = ImageFont.truetype(
                        '/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf',
                        font_size_px)
                except (IOError, OSError):
                    try:
                        font = ImageFont.truetype('DejaVuSans.ttf',
                                                  font_size_px)
                    except (IOError, OSError):
                        font = ImageFont.load_default()
                draw.text((x_px, y_px), data['text'], fill='black', font=font)

        pdf_surface.finish()
        print(f"Saved PDF (hybrid vector+raster): {pdf_path}")

        png_canvas.save(str(png_path), dpi=(self.RENDER_DPI, self.RENDER_DPI))
        print(f"Saved PNG preview: {png_path}")


def fractional_to_mm(frac_rect, fig_width_mm, fig_height_mm):
    """Convert matplotlib fractional [left, bottom, width, height] to
    (x_mm, y_mm_from_top, w_mm, h_mm).

    In matplotlib, (0,0) is bottom-left and fractions go up.
    In VectorAssembler, (0,0) is top-left and y increases downward.

    Args:
        frac_rect: [left, bottom, width, height] in fractional coords.
        fig_width_mm: Total figure width in mm.
        fig_height_mm: Total figure height in mm.

    Returns:
        (x_mm, y_mm, w_mm, h_mm) in top-left coordinate system.
    """
    left, bottom, width, height = frac_rect
    x_mm = left * fig_width_mm
    w_mm = width * fig_width_mm
    h_mm = height * fig_height_mm
    # Convert bottom-up to top-down: y_top = total - (bottom + height)
    y_mm = fig_height_mm - (bottom + height) * fig_height_mm
    return x_mm, y_mm, w_mm, h_mm
