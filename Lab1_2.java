import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;

public class Lab1_2 {

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

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.print("Enter input: ");
        String s = scanner.nextLine();
        System.out.println(percentageOf(s));
        scanner.close();
    }
}
