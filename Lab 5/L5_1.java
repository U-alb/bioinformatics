import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;

public class L5_1 {
    public static void main(String[] args) {
        String fastaFilePath = "rna.fna";
        try {
            String originalSequence = readFASTASequence(fastaFilePath);
            List<String> samples = generateRandomSamples(originalSequence);
            String reconstructedSequence = reconstructSequence(samples);
            writeAnalysisToFile(originalSequence, reconstructedSequence);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    static String readFASTASequence(String filePath) throws IOException {
        StringBuilder sequence = new StringBuilder();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    continue;
                }

                line = line.trim().toUpperCase();

                // Ensure only valid DNA nucleotides are included
                line = line.replaceAll("[^ATCG]", "");

                if (!line.isEmpty()) {
                    sequence.append(line);
                }
            }
        }

        return sequence.toString();
    }

    static List<String> generateRandomSamples(String sequence) {
        List<String> samples = new ArrayList<>();
        SecureRandom random = new SecureRandom();

        for (int i = 0; i < 2000; i++) {
            int sampleLength = 100 + random.nextInt(51);
            int startIndex = random.nextInt(sequence.length() - sampleLength);

            String sample = sequence.substring(startIndex, startIndex + sampleLength);
            samples.add(sample);
        }

        return samples;
    }

    private static String reconstructSequence(List<String> samples) {
        StringBuilder reconstructed = new StringBuilder();

        for (String sample : samples) {
            if (reconstructed.isEmpty()) {
                reconstructed.append(sample);
            } else {
                // Naive overlap detection
                for (int overlapLength = Math.min(sample.length(), reconstructed.length()); overlapLength > 0; overlapLength--) {
                    if (reconstructed.toString().endsWith(sample.substring(0, overlapLength))) {
                        reconstructed.append(sample.substring(overlapLength));
                        break;
                    }
                }
            }
        }

        return reconstructed.toString();
    }

    private static void writeAnalysisToFile(String original, String reconstructed) throws IOException {
        try (FileWriter writer = new FileWriter("answer.txt")) {
            writer.write(String.format("""
        DNA Sequence Reconstruction Analysis

        Original Sequence Length: %d nucleotides
        Reconstructed Sequence Length: %d nucleotides

        Sequence Reconstruction Similarity: %.2f%%

        Main Problems with the Reconstruction Approach:
        1. Incomplete Overlap Resolution:
           - Random sampling makes precise sequence reconstruction extremely challenging.
           - Overlapping regions may not capture full sequence context.

        2. Computational Complexity: Naive reconstruction method has exponential time complexity

        3. Sampling Bias:
           - Random sampling may not uniformly represent the entire sequence.
           - Some regions might be over- or under-represented.

        4. Short Read Limitations:
           - 100-150 base samples are too short for reliable reconstruction.
           - Repetitive sequences and structural variations complicate reassembly.

        5. Lack of Error Correction: No mechanism to handle potential sequencing errors or mutations.
        """,
                    original.length(),
                    reconstructed.length(),
                    calculateSequenceSimilarity(original, reconstructed) * 100)
            );
        }
    }

    // Calculate sequence similarity using Levenshtein distance
    private static double calculateSequenceSimilarity(String original, String reconstructed) {
        int maxLength = Math.max(original.length(), reconstructed.length());
        int editDistance = levenshteinDistance(original, reconstructed);
        return 1.0 - ((double) editDistance / maxLength);
    }
    
    private static int levenshteinDistance(String s1, String s2) {
        int[][] dp = new int[s1.length() + 1][s2.length() + 1];

        for (int i = 0; i <= s1.length(); i++) {
            for (int j = 0; j <= s2.length(); j++) {
                if (i == 0) {
                    dp[i][j] = j;
                } else if (j == 0) {
                    dp[i][j] = i;
                } else {
                    dp[i][j] = Math.min(
                            dp[i - 1][j] + 1,  // Deletion
                            Math.min(
                                    dp[i][j - 1] + 1,  // Insertion
                                    dp[i - 1][j - 1] + (s1.charAt(i - 1) == s2.charAt(j - 1) ? 0 : 1)  // Substitution
                            )
                    );
                }
            }
        }

        return dp[s1.length()][s2.length()];
    }
}
