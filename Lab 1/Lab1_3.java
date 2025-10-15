import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class Lab1_3 {

    public static String percentageOf(String s) {
        if (s == null || s.isEmpty()) {
            return "";
        }

        Map<Character, Integer> counts = new TreeMap<>();
        int total = 0;
        for (char c : s.toCharArray()) {
            counts.put(c, counts.getOrDefault(c, 0) + 1);
            total++;
        }

        TreeSet<Character> alphabet = new TreeSet<>(counts.keySet());
        StringBuilder builder = new StringBuilder();
        for (char c : alphabet) {
            builder.append(c);
        }

        StringBuilder percentBuilder = new StringBuilder();
        percentBuilder.append("Alphabet: ").append(builder).append(System.lineSeparator());
        percentBuilder.append("Percentages:").append(System.lineSeparator());
        for (char c : alphabet) {
            int cnt = counts.getOrDefault(c, 0);
            double pct = (cnt * 100.0) / total;
            percentBuilder.append(c)
                    .append(": ")
                    .append(String.format("%.2f", pct))
                    .append("%")
                    .append(System.lineSeparator());
        }

        return percentBuilder.toString().trim();
    }

    private static String readFastaSequences(File file) throws IOException {
        StringBuilder seqBuilder = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.startsWith(">")) continue; // skip header lines
                // append sequence lines (remove whitespace)
                seqBuilder.append(line.replaceAll("\\s+", ""));
            }
        }
        return seqBuilder.toString();
    }

    private static void createAndShowGui() {
        JFrame frame = new JFrame("FASTA Percent Composition");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(600, 400);
        frame.setLocationRelativeTo(null);

        JTextArea outputArea = new JTextArea();
        outputArea.setEditable(false);
        outputArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        JScrollPane scrollPane = new JScrollPane(outputArea);

        JButton openButton = new JButton("Open FASTA File");
        openButton.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setFileFilter(new FileNameExtensionFilter("FASTA files (*.fa, *.fasta, *.faa)", "fa", "fasta", "faa"));
            int ret = chooser.showOpenDialog(frame);
            if (ret == JFileChooser.APPROVE_OPTION) {
                File f = chooser.getSelectedFile();
                try {
                    String seq = readFastaSequences(f);
                    if (seq.isEmpty()) {
                        outputArea.setText("No sequence found in file.");
                    } else {
                        String result = percentageOf(seq);
                        outputArea.setText("File: " + f.getName() + System.lineSeparator() + System.lineSeparator() + result);
                    }
                } catch (IOException ex) {
                    outputArea.setText("Error reading file: " + ex.getMessage());
                }
            }
        });

        JPanel topPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        topPanel.add(openButton);

        frame.getContentPane().add(topPanel, BorderLayout.NORTH);
        frame.getContentPane().add(scrollPane, BorderLayout.CENTER);

        frame.setVisible(true);
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(Lab1_3::createAndShowGui);
    }
}
