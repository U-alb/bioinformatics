import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.List;

public class Lab2_3 {
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> new MainFrame().setVisible(true));
    }
}

class MainFrame extends JFrame {
    private ChartPanel chartPanel;
    private JLabel statusLabel;
    private static final int WINDOW_SIZE = 30;

    MainFrame() {
        super("Sliding Window Nucleotide Frequencies");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(900, 600);
        setLocationRelativeTo(null);

        JButton openBtn = new JButton("Open FASTA...");
        openBtn.addActionListener(e -> openFasta());

        statusLabel = new JLabel("Ready");

        chartPanel = new ChartPanel();

        JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT));
        top.add(openBtn);
        top.add(statusLabel);

        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(top, BorderLayout.NORTH);
        getContentPane().add(chartPanel, BorderLayout.CENTER);
    }

    private void openFasta() {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new FileNameExtensionFilter("FASTA files", "fa", "fasta", "fna", "ffn"));
        int rv = chooser.showOpenDialog(this);
        if (rv != JFileChooser.APPROVE_OPTION) return;
        File f = chooser.getSelectedFile();
        try {
            String seq = FastaParser.parse(Paths.get(f.toURI()));
            seq = seq.replaceAll("\\s+", "").toUpperCase();
            if (seq.length() < WINDOW_SIZE) {
                JOptionPane.showMessageDialog(this, "Sequence shorter than window size (" + WINDOW_SIZE + ").", "Too short", JOptionPane.INFORMATION_MESSAGE);
                return;
            }
            statusLabel.setText("Sequence length: " + seq.length());
            SlidingWindowAnalyzer.Result res = SlidingWindowAnalyzer.analyze(seq, WINDOW_SIZE);
            chartPanel.setData(res);
            chartPanel.repaint();
        } catch (IOException ex) {
            JOptionPane.showMessageDialog(this, "Error reading file: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
        }
    }
}

class FastaParser {
    // Reads all sequences in FASTA and concatenates them; ignores headers
    static String parse(Path path) throws IOException {
        StringBuilder sb = new StringBuilder();
        try (BufferedReader r = Files.newBufferedReader(path)) {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith(">")) continue;
                sb.append(line.trim());
            }
        }
        return sb.toString();
    }
}

class SlidingWindowAnalyzer {
    // Result holds alphabet in order and frequency vectors (one double[] per symbol)
    static class Result {
        List<String> alphabet; // e.g., ["A","C","G","T"]
        Map<String, double[]> freqMap; // symbol -> vector length = number of windows
        int windowSize;
    }

    static Result analyze(String seq, int windowSize) {
        int len = seq.length();
        int windows = Math.max(0, len - windowSize + 1);
        // Determine alphabet present (A,C,G,T prioritized)
        LinkedHashSet<Character> present = new LinkedHashSet<>();
        for (char c : seq.toCharArray()) {
            if (!Character.isWhitespace(c)) present.add(c);
        }
        // Build ordered alphabet: A,C,G,T first if present, then others in lexical order
        List<String> alpha = new ArrayList<>();
        for (char c : new char[] {'A','C','G','T'}) {
            if (present.contains(c)) alpha.add(String.valueOf(c));
        }
        for (char c : present) {
            String s = String.valueOf(c);
            if (!alpha.contains(s)) alpha.add(s);
        }
        // Limit to 4 signals max (per requirement); if fewer than 4, plot what's available
        if (alpha.size() > 4) alpha = alpha.subList(0, 4);

        Map<String, double[]> freqMap = new LinkedHashMap<>();
        for (String s : alpha) freqMap.put(s, new double[windows]);

        // Slide
        for (int i = 0; i < windows; i++) {
            int start = i;
            int end = i + windowSize; // exclusive
            // count in window
            Map<String, Integer> counts = new HashMap<>();
            for (int j = start; j < end; j++) {
                String k = String.valueOf(seq.charAt(j));
                counts.put(k, counts.getOrDefault(k, 0) + 1);
            }
            for (String sym : alpha) {
                int cnt = counts.getOrDefault(sym, 0);
                freqMap.get(sym)[i] = (double) cnt / windowSize;
            }
        }

        Result r = new Result();
        r.alphabet = alpha;
        r.freqMap = freqMap;
        r.windowSize = windowSize;
        return r;
    }
}

class ChartPanel extends JPanel {
    private SlidingWindowAnalyzer.Result data;
    private static final Color[] COLORS = { Color.RED, Color.BLUE, Color.GREEN.darker(), Color.MAGENTA };

    void setData(SlidingWindowAnalyzer.Result d) {
        this.data = d;
    }

    @Override
    protected void paintComponent(Graphics g0) {
        super.paintComponent(g0);
        Graphics2D g = (Graphics2D) g0;
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        if (data == null) {
            g.drawString("Open a FASTA file to display sliding-window frequencies.", 20, 20);
            return;
        }
        int w = getWidth();
        int h = getHeight();
        int left = 60, right = 20, top = 20, bottom = 60;
        int plotW = w - left - right;
        int plotH = h - top - bottom;
        // axes
        g.setColor(Color.WHITE);
        g.fillRect(left, top, plotW, plotH);
        g.setColor(Color.BLACK);
        g.drawRect(left, top, plotW, plotH);
        // labels
        g.drawString("Relative frequency", 8, top + 12);
        g.drawString("Windows", left + plotW/2 - 20, h - 20);

        int windows = data.freqMap.values().iterator().next().length;
        // y ticks 0.0..1.0
        for (int t = 0; t <= 10; t++) {
            int yy = top + plotH - (t * plotH / 10);
            g.setColor(new Color(220,220,220));
            g.drawLine(left, yy, left + plotW, yy);
            g.setColor(Color.BLACK);
            String lab = String.format("%.1f", t / 10.0);
            g.drawString(lab, 8, yy + 4);
        }

        // draw each signal
        int idx = 0;
        for (Map.Entry<String, double[]> e : data.freqMap.entrySet()) {
            double[] vec = e.getValue();
            g.setColor(COLORS[idx % COLORS.length]);
            // label legend
            g.fillRect(w - 120, top + 10 + idx*18, 12, 12);
            g.setColor(Color.BLACK);
            g.drawString(e.getKey(), w - 100, top + 10 + idx*18 + 12);
            // draw polyline
            g.setColor(COLORS[idx % COLORS.length]);
            int prevX = -1, prevY = -1;
            for (int i = 0; i < vec.length; i++) {
                double v = vec[i]; // 0..1
                int x = left + (int) Math.round((double)i * (plotW - 1) / Math.max(1, windows - 1));
                int y = top + plotH - (int) Math.round(v * plotH);
                if (prevX != -1) {
                    g.drawLine(prevX, prevY, x, y);
                }
                prevX = x; prevY = y;
            }
            idx++;
        }
    }
}
