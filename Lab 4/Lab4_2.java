import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

import javax.swing.*;
import java.io.*;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Lab4_2 {
    private static final Map<String, String> geneticCode = new HashMap<>();

    static {
        geneticCode.put("UUU", "Phe");
        geneticCode.put("UUC", "Phe");
        geneticCode.put("UUA", "Leu");
        geneticCode.put("UUG", "Leu");
        geneticCode.put("UCU", "Ser");
        geneticCode.put("UCC", "Ser");
        geneticCode.put("UCA", "Ser");
        geneticCode.put("UCG", "Ser");
        geneticCode.put("UAU", "Tyr");
        geneticCode.put("UAC", "Tyr");
        geneticCode.put("UAA", "Stop");
        geneticCode.put("UAG", "Stop");
        geneticCode.put("UGU", "Cys");
        geneticCode.put("UGC", "Cys");
        geneticCode.put("UGA", "Stop");
        geneticCode.put("UGG", "Trp");
        geneticCode.put("CUU", "Leu");
        geneticCode.put("CUC", "Leu");
        geneticCode.put("CUA", "Leu");
        geneticCode.put("CUG", "Leu");
        geneticCode.put("CCU", "Pro");
        geneticCode.put("CCC", "Pro");
        geneticCode.put("CCA", "Pro");
        geneticCode.put("CCG", "Pro");
        geneticCode.put("CAU", "His");
        geneticCode.put("CAC", "His");
        geneticCode.put("CAA", "Gln");
        geneticCode.put("CAG", "Gln");
        geneticCode.put("CGU", "Arg");
        geneticCode.put("CGC", "Arg");
        geneticCode.put("CGA", "Arg");
        geneticCode.put("CGG", "Arg");
        geneticCode.put("AUU", "Ile");
        geneticCode.put("AUC", "Ile");
        geneticCode.put("AUA", "Ile");
        geneticCode.put("AUG", "Met");
        geneticCode.put("ACU", "Thr");
        geneticCode.put("ACC", "Thr");
        geneticCode.put("ACA", "Thr");
        geneticCode.put("ACG", "Thr");
        geneticCode.put("AAU", "Asn");
        geneticCode.put("AAC", "Asn");
        geneticCode.put("AAA", "Lys");
        geneticCode.put("AAG", "Lys");
        geneticCode.put("AGU", "Ser");
        geneticCode.put("AGC", "Ser");
        geneticCode.put("AGA", "Arg");
        geneticCode.put("AGG", "Arg");
        geneticCode.put("GUU", "Val");
        geneticCode.put("GUC", "Val");
        geneticCode.put("GUA", "Val");
        geneticCode.put("GUG", "Val");
        geneticCode.put("GCU", "Ala");
        geneticCode.put("GCC", "Ala");
        geneticCode.put("GCA", "Ala");
        geneticCode.put("GCG", "Ala");
        geneticCode.put("GAU", "Asp");
        geneticCode.put("GAC", "Asp");
        geneticCode.put("GAA", "Glu");
        geneticCode.put("GAG", "Glu");
        geneticCode.put("GGU", "Gly");
        geneticCode.put("GGC", "Gly");
        geneticCode.put("GGA", "Gly");
        geneticCode.put("GGG", "Gly");
    }

    public static void main(String[] args) {
        try {
            // Read COVID-19 and Influenza genome sequences
            String covidSequence = readFastaFile("covid19.fasta");
            String influenzaSequence = readFastaFile("influenza-river-sequence.fasta");

            // Analyze and create charts for COVID-19 and Influenza genomes
            createCodonFrequencyChart(covidSequence, "COVID-19 Top 10 Codons");
            createCodonFrequencyChart(influenzaSequence, "Influenza Top 10 Codons");

            // Task c: Compare codons between the two genomes
            compareGenomeCodons(covidSequence, influenzaSequence);

            // Task d: Show top 3 amino acids for each genome
            analyzeAminoAcidFrequencies(covidSequence, "COVID-19");
            analyzeAminoAcidFrequencies(influenzaSequence, "Influenza");

            // Additional analysis
            performDetailedGenomeAnalysis(covidSequence, influenzaSequence);

            // Codon usage bias comparison
            compareCodonUsageBias(covidSequence, influenzaSequence);

        } catch (IOException e) {
            System.err.println("Error reading genome files: " + e.getMessage());
        }
    }
    private static void createCodonFrequencyChart(String sequence, String chartTitle) {
        // Transcribe DNA to RNA
        String rnaSequence = sequence.replace('T', 'U');

        // Count codon frequencies
        Map<String, Integer> codonFrequencies = new HashMap<>();
        for (int i = 0; i < rnaSequence.length() - 2; i += 3) {
            String codon = rnaSequence.substring(i, i + 3);
            codonFrequencies.put(codon, codonFrequencies.getOrDefault(codon, 0) + 1);
        }

        // Get top 10 most frequent codons
        List<Map.Entry<String, Integer>> topCodons = codonFrequencies.entrySet().stream()
                .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                .limit(10)
                .toList();

        // Create dataset for the chart
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (Map.Entry<String, Integer> entry : topCodons) {
            dataset.addValue(entry.getValue(), "Frequency", entry.getKey());
        }

        // Create the chart
        JFreeChart chart = ChartFactory.createBarChart(
                chartTitle,       // Chart title
                "Codons",         // X-Axis Label
                "Frequency",      // Y-Axis Label
                dataset,          // Dataset
                PlotOrientation.VERTICAL, // Plot Orientation
                false,            // Show Legend
                true,             // Use tooltips
                false             // Configure chart to generate URLs?
        );

        // Create and display the chart
        ChartPanel chartPanel = new ChartPanel(chart);
        JFrame frame = new JFrame(chartTitle);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.setSize(800, 600);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    private static String readFastaFile(String filename) throws IOException {
        // Print out the full path for debugging
        System.out.println("Attempting to read file: " + filename);

        // Get the input stream and print diagnostic information
        InputStream inputStream = Lab4_2.class.getClassLoader().getResourceAsStream(filename);

        if (inputStream == null) {
            // Detailed error logging
            System.err.println("File not found in resources: " + filename);
            System.err.println("Current working directory: " + System.getProperty("user.dir"));
            System.err.println("Classpath: " + System.getProperty("java.class.path"));

            // List all resources in the classpath
            try {
                URL[] urls = ((URLClassLoader)Lab4_2.class.getClassLoader()).getURLs();
                for (URL url : urls) {
                    System.err.println("Classpath entry: " + url);
                }
            } catch (Exception e) {
                System.err.println("Could not list classpath entries");
            }

            throw new FileNotFoundException("Resource file not found: " + filename);
        }

        StringBuilder sequence = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(inputStream))) {
            String line;
            boolean isSequenceLine = false;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    isSequenceLine = true;
                    continue;
                }
                if (isSequenceLine) {
                    sequence.append(line.trim().toUpperCase());
                }
            }
        }
        return sequence.toString();
    }


    // Optional: Method to print top 10 codons to console
    private static void printTopCodons(Map<String, Integer> codonFrequencies, String genomeName) {
        System.out.println("\nTop 10 Most Frequent Codons for " + genomeName + ":");
        codonFrequencies.entrySet().stream()
                .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                .limit(10)
                .forEach(entry -> System.out.println(entry.getKey() + ": " + entry.getValue()));
    }

    // Additional method to compare codons between genomes
    private static void compareGenomeCodons(String covidSequence, String influenzaSequence) {
        // Transcribe both sequences to RNA
        String covidRnaSequence = covidSequence.replace('T', 'U');
        String influenzaRnaSequence = influenzaSequence.replace('T', 'U');

        // Count codon frequencies for both genomes
        Map<String, Integer> covidCodonFrequencies = getCodonFrequencies(covidRnaSequence);
        Map<String, Integer> influenzaCodonFrequencies = getCodonFrequencies(influenzaRnaSequence);

        // Create a comparison chart
        createComparativeCodonChart(covidCodonFrequencies, influenzaCodonFrequencies);
    }

    private static void createComparativeCodonChart(
            Map<String, Integer> covidCodonFrequencies,
            Map<String, Integer> influenzaCodonFrequencies) {

        // Combine and get top 10 unique codons across both genomes
        Set<String> topCodons = Stream.concat(
                        covidCodonFrequencies.entrySet().stream(),
                        influenzaCodonFrequencies.entrySet().stream())
                .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                .limit(10)
                .map(Map.Entry::getKey)
                .collect(Collectors.toSet());

        // Create dataset for the comparative chart
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (String codon : topCodons) {
            int covidFreq = covidCodonFrequencies.getOrDefault(codon, 0);
            int influenzaFreq = influenzaCodonFrequencies.getOrDefault(codon, 0);

            dataset.addValue(covidFreq, "COVID-19", codon);
            dataset.addValue(influenzaFreq, "Influenza", codon);
        }

        // Create the comparative chart
        JFreeChart chart = ChartFactory.createBarChart(
                "Comparative Codon Frequencies",  // Chart title
                "Codons",                         // X-Axis Label
                "Frequency",                      // Y-Axis Label
                dataset,                          // Dataset
                PlotOrientation.VERTICAL,         // Plot Orientation
                true,                             // Show Legend
                true,                             // Use tooltips
                false                             // Configure chart to generate URLs?
        );

        // Create and display the chart
        ChartPanel chartPanel = new ChartPanel(chart);
        JFrame frame = new JFrame("Comparative Codon Frequencies");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.setSize(1000, 600);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    private static Map<String, Integer> getCodonFrequencies(String rnaSequence) {
        Map<String, Integer> codonFrequencies = new HashMap<>();
        for (int i = 0; i < rnaSequence.length() - 2; i += 3) {
            String codon = rnaSequence.substring(i, i + 3);
            codonFrequencies.put(codon, codonFrequencies.getOrDefault(codon, 0) + 1);
        }
        return codonFrequencies;
    }

    // Method to analyze amino acid frequencies
    private static void analyzeAminoAcidFrequencies(String sequence, String genomeName) {
        // Transcribe DNA to RNA
        String rnaSequence = sequence.replace('T', 'U');

        // Translate to amino acids and count frequencies
        Map<String, Integer> aminoAcidFrequencies = new HashMap<>();
        for (int i = 0; i < rnaSequence.length() - 2; i += 3) {
            String codon = rnaSequence.substring(i, i + 3);
            String aminoAcid = geneticCode.get(codon);
            if (aminoAcid != null && !aminoAcid.equals("Stop")) {
                aminoAcidFrequencies.put(aminoAcid,
                        aminoAcidFrequencies.getOrDefault(aminoAcid, 0) + 1);
            }
        }

        // Create amino acid frequency chart
        createAminoAcidFrequencyChart(aminoAcidFrequencies, genomeName);
    }

    private static void createAminoAcidFrequencyChart(
            Map<String, Integer> aminoAcidFrequencies,
            String genomeName) {

        // Get top 3 amino acids
        List<Map.Entry<String, Integer>> topAminoAcids = aminoAcidFrequencies.entrySet().stream()
                .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                .limit(3)
                .toList();

        // Create dataset for the chart
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (Map.Entry<String, Integer> entry : topAminoAcids) {
            dataset.addValue(entry.getValue(), "Frequency", entry.getKey());
        }

        // Create the chart
        JFreeChart chart = ChartFactory.createBarChart(
                genomeName + " - Top 3 Amino Acids",  // Chart title
                "Amino Acids",                        // X-Axis Label
                "Frequency",                          // Y-Axis Label
                dataset,                              // Dataset
                PlotOrientation.VERTICAL,             // Plot Orientation
                false,                                // Show Legend
                true,                                 // Use tooltips
                false                                 // Configure chart to generate URLs?
        );

        // Create and display the chart
        ChartPanel chartPanel = new ChartPanel(chart);
        JFrame frame = new JFrame(genomeName + " Amino Acid Frequencies");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.setSize(800, 600);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    // Detailed genome analysis method
    private static void performDetailedGenomeAnalysis(
            String covidSequence,
            String influenzaSequence) {

        // Analyze codon frequencies
        Map<String, Integer> covidCodonFreq = getCodonFrequencies(
                covidSequence.replace('T', 'U')
        );
        Map<String, Integer> influenzaCodonFreq = getCodonFrequencies(
                influenzaSequence.replace('T', 'U')
        );

        // Print top codons to console
        System.out.println("COVID-19 Top Codons:");
        printTopCodons(covidCodonFreq, "COVID-19");

        System.out.println("\nInfluenza Top Codons:");
        printTopCodons(influenzaCodonFreq, "Influenza");

        // Create visualizations
        createCodonFrequencyChart(covidSequence, "COVID-19 Top 10 Codons");
        createCodonFrequencyChart(influenzaSequence, "Influenza Top 10 Codons");

        // Create comparative codon chart
        createComparativeCodonChart(covidCodonFreq, influenzaCodonFreq);

        // Analyze and visualize amino acid frequencies
        analyzeAminoAcidFrequencies(covidSequence, "COVID-19");
        analyzeAminoAcidFrequencies(influenzaSequence, "Influenza");
    }

    // Additional utility methods for advanced analysis
    private static Map<String, Double> calculateCodonUsageBias(String sequence) {
        String rnaSequence = sequence.replace('T', 'U');
        Map<String, Integer> codonFrequencies = getCodonFrequencies(rnaSequence);
        int totalCodons = codonFrequencies.values().stream().mapToInt(Integer::intValue).sum();

        return codonFrequencies.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        entry -> (double) entry.getValue() / totalCodons
                ));
    }

    private static void compareCodonUsageBias(String covidSequence, String influenzaSequence) {
        Map<String, Double> covidCodonBias = calculateCodonUsageBias(covidSequence);
        Map<String, Double> influenzaCodonBias = calculateCodonUsageBias(influenzaSequence);

        // Create a dataset for comparing codon usage bias
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        // Combine unique codons from both sequences
        Set<String> uniqueCodons = new HashSet<>(covidCodonBias.keySet());
        uniqueCodons.addAll(influenzaCodonBias.keySet());

        // Add data to dataset
        for (String codon : uniqueCodons) {
            dataset.addValue(
                    covidCodonBias.getOrDefault(codon, 0.0),
                    "COVID-19",
                    codon
            );
            dataset.addValue(
                    influenzaCodonBias.getOrDefault(codon, 0.0),
                    "Influenza",
                    codon
            );
        }

        // Create and display the chart
        JFreeChart chart = ChartFactory.createBarChart(
                "Codon Usage Bias Comparison",  // Chart title
                "Codons",                       // X-Axis Label
                "Usage Frequency",              // Y-Axis Label
                dataset,                        // Dataset
                PlotOrientation.VERTICAL,       // Plot Orientation
                true,                           // Show Legend
                true,                           // Use tooltips
                false                           // Configure chart to generate URLs?
        );

        // Display chart
        ChartPanel chartPanel = new ChartPanel(chart);
        JFrame frame = new JFrame("Codon Usage Bias Comparison");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.setSize(1200, 800);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}
