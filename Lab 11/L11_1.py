import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox

def needleman_wunsch_matrix(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    ptr = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
        ptr[i][0] = 'U'
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap
        ptr[0][j] = 'L'

    for i in range(1, n + 1):
        a = seq1[i - 1]
        for j in range(1, m + 1):
            b = seq2[j - 1]
            if a == b:
                diag = score[i - 1][j - 1] + match
            else:
                diag = score[i - 1][j - 1] + mismatch
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diag, up, left)
            score[i][j] = best
            if best == diag:
                ptr[i][j] = 'D'
            elif best == up:
                ptr[i][j] = 'U'
            else:
                ptr[i][j] = 'L'

    return score, ptr

def traceback_alignment(seq1, seq2, ptr):
    i = len(seq1)
    j = len(seq2)
    path = [(i, j)]
    while i > 0 or j > 0:
        p = ptr[i][j]
        if p == 'D':
            i -= 1
            j -= 1
        elif p == 'U':
            i -= 1
        elif p == 'L':
            j -= 1
        else:
            if i > 0:
                i -= 1
            else:
                j -= 1
        path.append((i, j))
    path.reverse()
    return path

def build_alignment_from_path(seq1, seq2, path):
    align1 = []
    align2 = []
    matches = 0
    for k in range(1, len(path)):
        i0, j0 = path[k - 1]
        i1, j1 = path[k]
        if i1 == i0 + 1 and j1 == j0 + 1:
            a = seq1[i0]
            b = seq2[j0]
            align1.append(a)
            align2.append(b)
            if a == b:
                matches += 1
        elif i1 == i0 + 1 and j1 == j0:
            align1.append(seq1[i0])
            align2.append('-')
        elif i1 == i0 and j1 == j0 + 1:
            align1.append('-')
            align2.append(seq2[j0])
        else:
            pass
    return ''.join(align1), ''.join(align2), matches

def hex_interp(c1, c2, t):
    def to_int(h): return int(h, 16)
    r1 = to_int(c1[1:3]); g1 = to_int(c1[3:5]); b1 = to_int(c1[5:7])
    r2 = to_int(c2[1:3]); g2 = to_int(c2[3:5]); b2 = to_int(c2[5:7])
    r = int(r1 + (r2 - r1) * t)
    g = int(g1 + (g2 - g1) * t)
    b = int(b1 + (b2 - b1) * t)
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)

def score_to_color(val, vmin, vmax):
    if vmax == vmin:
        t = 0.5
    else:
        t = (val - vmin) / (vmax - vmin)
        t = max(0.0, min(1.0, t))
    if t < 0.5:
        return hex_interp('#17202a', '#9b2b6f', t * 2)
    else:
        return hex_interp('#9b2b6f', '#ff3b30', (t - 0.5) * 2)

class NWGui:
    def __init__(self, root):
        self.root = root
        root.title("Needleman–Wunsch Alignment (Tkinter)")
        self._build_ui()
        # default sequences
        self.seq1_entry.insert(0, "ACCGTGAAGCCAATAC")
        self.seq2_entry.insert(0, "AGCGTGCAGCCAATAC")

    def _build_ui(self):
        frm = ttk.Frame(self.root, padding=8)
        frm.grid(sticky='nsew')
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        left = ttk.Frame(frm)
        left.grid(row=0, column=0, sticky='nw', padx=(0,10))

        ttk.Label(left, text="Seq 1:").grid(row=0, column=0, sticky='w')
        self.seq1_entry = ttk.Entry(left, width=40)
        self.seq1_entry.grid(row=0, column=1, sticky='w')

        ttk.Label(left, text="Seq 2:").grid(row=1, column=0, sticky='w')
        self.seq2_entry = ttk.Entry(left, width=40)
        self.seq2_entry.grid(row=1, column=1, sticky='w')

        score_frame = ttk.LabelFrame(left, text="Scoring")
        score_frame.grid(row=2, column=0, columnspan=2, pady=(8,8), sticky='w')
        ttk.Label(score_frame, text="Match:").grid(row=0, column=0)
        self.match_spin = tk.Spinbox(score_frame, from_=-5, to=10, width=5)
        self.match_spin.delete(0,"end"); self.match_spin.insert(0,"1")
        self.match_spin.grid(row=0, column=1, padx=4)
        ttk.Label(score_frame, text="Mismatch:").grid(row=0, column=2)
        self.mismatch_spin = tk.Spinbox(score_frame, from_=-10, to=5, width=5)
        self.mismatch_spin.delete(0,"end"); self.mismatch_spin.insert(0,"-1")
        self.mismatch_spin.grid(row=0, column=3, padx=4)
        ttk.Label(score_frame, text="Gap:").grid(row=0, column=4)
        self.gap_spin = tk.Spinbox(score_frame, from_=-10, to=0, width=5)
        self.gap_spin.delete(0,"end"); self.gap_spin.insert(0,"-2")
        self.gap_spin.grid(row=0, column=5, padx=4)

        opt_frame = ttk.Frame(left)
        opt_frame.grid(row=3, column=0, columnspan=2, pady=(4,8), sticky='w')
        self.show_numbers_var = tk.BooleanVar(value=True)
        chk = ttk.Checkbutton(opt_frame, text="Show numbers in cells", variable=self.show_numbers_var)
        chk.grid(row=0, column=0, sticky='w')

        self.align_btn = ttk.Button(left, text="Align", command=self.run_alignment)
        self.align_btn.grid(row=4, column=0, columnspan=2, pady=(4,4))

        canv_frame = ttk.Frame(frm)
        canv_frame.grid(row=0, column=1, sticky='n', padx=(10,0))

        self.matrix_canvas = tk.Canvas(canv_frame, width=420, height=420, bg='white')
        self.matrix_canvas.grid(row=0, column=0, padx=(0,8), pady=(0,8))
        self.trace_canvas = tk.Canvas(canv_frame, width=420, height=420, bg='white')
        self.trace_canvas.grid(row=0, column=1, pady=(0,8))

        bottom = ttk.Frame(frm)
        bottom.grid(row=1, column=0, columnspan=2, sticky='nsew')
        ttk.Label(bottom, text="Alignment result:").grid(row=0, column=0, sticky='w')
        self.result_txt = scrolledtext.ScrolledText(bottom, height=8, width=100, wrap='none', font=("Courier", 10))
        self.result_txt.grid(row=1, column=0, pady=(4,4))

    def run_alignment(self):
        seq1 = self.seq1_entry.get().strip().upper()
        seq2 = self.seq2_entry.get().strip().upper()
        if not seq1 or not seq2:
            messagebox.showwarning("Input required", "Please enter both sequences.")
            return
        try:
            match = int(self.match_spin.get())
            mismatch = int(self.mismatch_spin.get())
            gap = int(self.gap_spin.get())
        except ValueError:
            messagebox.showerror("Invalid scores", "Match/Mismatch/Gap must be integers.")
            return

        score_mat, ptr = needleman_wunsch_matrix(seq1, seq2, match=match, mismatch=mismatch, gap=gap)
        path = traceback_alignment(seq1, seq2, ptr)
        al1, al2, matches = build_alignment_from_path(seq1, seq2, path)

        alignment_score = score_mat[len(seq1)][len(seq2)]
        sim_pct = 100.0 * matches / max(1, len(al1))

        self.draw_score_matrix(score_mat, seq1, seq2)
        self.draw_traceback_grid(path, seq1, seq2)

        self.result_txt.delete('1.0', tk.END)
        match_line = []
        for a, b in zip(al1, al2):
            if a == b and a != '-':
                match_line.append('|')
            elif a == '-' or b == '-':
                match_line.append(' ')
            else:
                match_line.append('.')
        match_line = ''.join(match_line)

        self.result_txt.insert(tk.END, f"Score: {alignment_score}\n")
        self.result_txt.insert(tk.END, f"Matches: {matches} / {len(al1)}  ({sim_pct:.1f}%)\n\n")
        self.result_txt.insert(tk.END, al1 + "\n")
        self.result_txt.insert(tk.END, match_line + "\n")
        self.result_txt.insert(tk.END, al2 + "\n")

    def draw_score_matrix(self, score_mat, seq1, seq2):
        c = self.matrix_canvas
        c.delete('all')
        n = len(seq1)
        m = len(seq2)
        rows = n + 1
        cols = m + 1
        max_dim = max(rows, cols)
        W = int(c['width'])
        H = int(c['height'])
        cell_size = max(12, min(30, min(W // cols, H // rows)))
        flat = [v for row in score_mat for v in row]
        vmin = min(flat)
        vmax = max(flat)

        show_numbers = self.show_numbers_var.get() and (rows * cols <= 800)
        font_size = 8 if cell_size >= 18 else 6
        for i in range(rows):
            for j in range(cols):
                x0 = j * cell_size
                y0 = i * cell_size
                x1 = x0 + cell_size
                y1 = y0 + cell_size
                val = score_mat[i][j]
                color = score_to_color(val, vmin, vmax)
                c.create_rectangle(x0, y0, x1, y1, fill=color, outline='gray')
                if show_numbers and cell_size >= 14:
                    c.create_text((x0 + x1) // 2, (y0 + y1) // 2,
                                  text=str(val), font=("Helvetica", font_size))

        for j in range(1, cols):
            x = j * cell_size + cell_size // 2
            y = cell_size // 8
            c.create_text(x, y, text=seq2[j - 1], anchor='n', font=("Helvetica", max(8, font_size)))
        for i in range(1, rows):
            x = cell_size // 8
            y = i * cell_size + cell_size // 2
            c.create_text(x, y, text=seq1[i - 1], anchor='w', font=("Helvetica", max(8, font_size)))

        for i in range(rows + 1):
            c.create_line(0, i * cell_size, cols * cell_size, i * cell_size, fill='#444444')
        for j in range(cols + 1):
            c.create_line(j * cell_size, 0, j * cell_size, rows * cell_size, fill='#444444')

    def draw_traceback_grid(self, path, seq1, seq2):
        c = self.trace_canvas
        c.delete('all')
        n = len(seq1)
        m = len(seq2)
        rows = n + 1
        cols = m + 1
        W = int(c['width'])
        H = int(c['height'])
        cell_size = max(12, min(30, min(W // cols, H // rows)))

        for i in range(rows):
            for j in range(cols):
                x0 = j * cell_size
                y0 = i * cell_size
                x1 = x0 + cell_size
                y1 = y0 + cell_size
                c.create_rectangle(x0, y0, x1, y1, fill='#fff5e6', outline='#cccccc')

        for (i, j) in path:
            x0 = j * cell_size
            y0 = i * cell_size
            x1 = x0 + cell_size
            y1 = y0 + cell_size
            c.create_rectangle(x0 + 1, y0 + 1, x1 - 1, y1 - 1, fill='#ff4b4b', outline='black')

        if path:
            s_i, s_j = path[0]
            e_i, e_j = path[-1]
            c.create_text((s_j + 0.5) * cell_size, (s_i + 0.5) * cell_size, text='S', font=("Helvetica", 10, "bold"))
            c.create_text((e_j + 0.5) * cell_size, (e_i + 0.5) * cell_size, text='E', font=("Helvetica", 10, "bold"))

        for j in range(1, cols):
            x = j * cell_size + cell_size // 2
            y = cell_size // 8
            c.create_text(x, y, text=seq2[j - 1], anchor='n', font=("Helvetica", 9))
        for i in range(1, rows):
            x = cell_size // 8
            y = i * cell_size + cell_size // 2
            c.create_text(x, y, text=seq1[i - 1], anchor='w', font=("Helvetica", 9))

        for i in range(rows + 1):
            c.create_line(0, i * cell_size, cols * cell_size, i * cell_size, fill='#444444')
        for j in range(cols + 1):
            c.create_line(j * cell_size, 0, j * cell_size, rows * cell_size, fill='#444444')


def main():
    root = tk.Tk()
    app = NWGui(root)
    root.mainloop()

if __name__ == "__main__":
    main()
