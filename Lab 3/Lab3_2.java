import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Lab3_2 extends Component {

    public static void main(String[] args) {
        // Launch the GUI
        SwingUtilities.invokeLater(() -> {
            new MeltingTemperatureGUI().setVisible(true);
        });
    }

    public static class MeltingTemperatureGUI extends JFrame {
        // GUI Components
        private JTextField filePathField;
        private JButton browseButton;
        private JButton calculateButton;
        private final JTable resultsTable;
        private JSpinner windowSizeSpinner;
        private final JFrame parentFrame = this;

        private JPanel createInputPanel() {
            JPanel inputPanel = new JPanel();
            inputPanel.setLayout(new FlowLayout());

            // File path input
            filePathField = new JTextField(30);
            filePathField.setEditable(false);
            inputPanel.add(new JLabel("FASTA File:"));
            inputPanel.add(filePathField);

            // Browse button
            browseButton = new JButton("Browse");
            browseButton.addActionListener(e -> browseFASTAFile());
            inputPanel.add(browseButton);

            // Window size spinner
            SpinnerNumberModel windowSizeModel = new SpinnerNumberModel(8, 1, 50, 1);
            windowSizeSpinner = new JSpinner(windowSizeModel);
            inputPanel.add(new JLabel("Window Size:"));
            inputPanel.add(windowSizeSpinner);

            // Calculate button
            calculateButton = new JButton("Calculate Melting Temperatures");
            calculateButton.addActionListener(e -> calculateMeltingTemperatures());
            calculateButton.setEnabled(false);
            inputPanel.add(calculateButton);

            return inputPanel;
        }

        private void browseFASTAFile() {
            JFileChooser fileChooser = new JFileChooser();
            fileChooser.setDialogTitle("Select FASTA File");

            // Set file filter for FASTA files
            FileNameExtensionFilter fastaFilter = new FileNameExtensionFilter(
                    "FASTA Files (*.fasta, *.fa)", "fasta", "fa"
            );
            fileChooser.setFileFilter(fastaFilter);

            int result = fileChooser.showOpenDialog(this);
            if (result == JFileChooser.APPROVE_OPTION) {
                File selectedFile = fileChooser.getSelectedFile();
                filePathField.setText(selectedFile.getAbsolutePath());
                calculateButton.setEnabled(true);
            }
        }

        private void calculateMeltingTemperatures() {
            // Add a progress indicator
            JDialog progressDialog = new JDialog(this, "Calculating", true);
            JProgressBar progressBar = new JProgressBar();
            progressBar.setIndeterminate(true);
            progressDialog.add(progressBar);
            progressDialog.setSize(300, 100);
            progressDialog.setLocationRelativeTo(this);

            // Run calculation in a background thread
            SwingWorker<List<MeltingTemperatureResult>, Void> worker = new SwingWorker<>() {
                @Override
                protected List<MeltingTemperatureResult> doInBackground() throws Exception {
                    try {
                        // Read the sequence from the FASTA file
                        String sequence = readFASTAFile(filePathField.getText());

                        // Get window size from spinner
                        int windowSize = (Integer) windowSizeSpinner.getValue();

                        // Calculate melting temperatures
                        return calculateSlidingWindowMeltingTemperatures(sequence, windowSize);

                    } catch (IOException ex) {
                        throw new Exception("Error reading FASTA file: " + ex.getMessage());
                    }
                }

                @Override
                protected void done() {
                    progressDialog.dispose();
                    try {
                        List<MeltingTemperatureResult> results = get();
                        // Display results in table
                        displayResults(results);
                    } catch (Exception ex) {
                        JOptionPane.showMessageDialog(
                                parentFrame,
                                ex.getMessage(),
                                "Error",
                                JOptionPane.ERROR_MESSAGE
                        );
                    }
                }
            };

            // Start the worker and show progress dialog
            worker.execute();
            progressDialog.setVisible(true);
        }

        public MeltingTemperatureGUI() {
            // Set up the main frame
            setTitle("DNA Melting Temperature Analysis");
            setSize(1000, 600);
            setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            setLayout(new BorderLayout());

            // Create input panel
            JPanel inputPanel = createInputPanel();
            add(inputPanel, BorderLayout.NORTH);

            // Create results table
            resultsTable = new JTable();
            JScrollPane scrollPane = new JScrollPane(resultsTable);
            add(scrollPane, BorderLayout.CENTER);

            // Center the frame
            setLocationRelativeTo(null);
        }

        private void displayResults(List<MeltingTemperatureResult> results) {
            // Create table model with column names
            String[] columnNames = {
                    "Start Position",
                    "Window Sequence",
                    "Simple Tm (°C)",
                    "Complex Tm (°C)"
            };

            // Convert results to 2D object array for table
            Object[][] data = new Object[results.size()][4];
            for (int i = 0; i < results.size(); i++) {
                MeltingTemperatureResult result = results.get(i);
                data[i] = new Object[]{
                        result.startPosition,
                        result.window,
                        String.format("%.2f", result.simpleMeltingTemperature),
                        String.format("%.2f", result.complexMeltingTemperature)
                };
            }

            // Create table model with non-editable cells
            DefaultTableModel tableModel = new DefaultTableModel(data, columnNames) {
                @Override
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };

            // Update results table
            resultsTable.setModel(tableModel);

            // Optional: Adjust column widths
            resultsTable.getColumnModel().getColumn(0).setPreferredWidth(50);
            resultsTable.getColumnModel().getColumn(1).setPreferredWidth(150);
            resultsTable.getColumnModel().getColumn(2).setPreferredWidth(50);
            resultsTable.getColumnModel().getColumn(3).setPreferredWidth(50);
        }
    }

    private static String readFASTAFile(String filePath) throws IOException {
        StringBuilder sequence = new StringBuilder();
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                // Skip header lines (starting with '>')
                if (!line.startsWith(">")) {
                    sequence.append(line.trim().toUpperCase());
                }
            }
        }
        return sequence.toString();
    }

    private static List<MeltingTemperatureResult> calculateSlidingWindowMeltingTemperatures(String sequence, int windowSize) {
        List<MeltingTemperatureResult> results = new ArrayList<>();

        // Ensure we don't go beyond the sequence length
        for (int i = 0; i <= sequence.length() - windowSize; i++) {
            // Extract window
            String window = sequence.substring(i, i + windowSize);

            // Count base frequencies
            int freqA = countBase(window, 'A');
            int freqT = countBase(window, 'T');
            int freqC = countBase(window, 'C');
            int freqG = countBase(window, 'G');

            // Calculate melting temperatures
            double simpleTm = calculateMeltingTemperatureSimple(freqA, freqT, freqC, freqG);
            double complexTm = calculateMeltingTemperatureComplex(window, freqA, freqT, freqC, freqG);

            // Store results
            results.add(new MeltingTemperatureResult(i, window, simpleTm, complexTm));
        }

        return results;
    }

    private static int countBase(String dna, char base) {
        return Lab3_1.countBase(dna, base);
    }

    private static double calculateMeltingTemperatureSimple(
            int freqA, int freqT, int freqC, int freqG
    ) {
        return 4 * (freqG + freqC) + 2 * (freqA + freqT);
    }

    private static double calculateMeltingTemperatureComplex(
            String dna, int freqA, int freqT, int freqC, int freqG
    ) {
        return Lab3_1.calculateMeltingTemperatureComplex(dna, 0, 0, freqC, freqG);
    }

    private static class MeltingTemperatureResult {
        int startPosition;
        String window;
        double simpleMeltingTemperature;
        double complexMeltingTemperature;

        public MeltingTemperatureResult(
                int startPosition,
                String window,
                double simpleMeltingTemperature,
                double complexMeltingTemperature
        ) {
            this.startPosition = startPosition;
            this.window = window;
            this.simpleMeltingTemperature = simpleMeltingTemperature;
            this.complexMeltingTemperature = complexMeltingTemperature;
        }
    }
}
