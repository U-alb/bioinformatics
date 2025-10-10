import java.util.*;

public class Lab2_2 {
    public static void main(String[] args) {
        String S = "ATTGTCCCAATCTGTTG";
        S = S.toUpperCase();

        Map<String,Integer> diCounts = new LinkedHashMap<>();
        Map<String,Integer> triCounts = new LinkedHashMap<>();

        int len = S.length();
        int diWindows = Math.max(0, len - 1);
        int triWindows = Math.max(0, len - 2);

        for (int i = 0; i < len; i++) {
            if (i + 2 <= len) {
                String di = S.substring(i, i + 2);
                diCounts.put(di, diCounts.getOrDefault(di, 0) + 1);
            }
            if (i + 3 <= len) {
                String tri = S.substring(i, i + 3);
                triCounts.put(tri, triCounts.getOrDefault(tri, 0) + 1);
            }
        }

        System.out.println("Sequence: " + S);
        System.out.println();

        System.out.println("Observed dinucleotides (count, frequency = count/" + diWindows + "):");
        for (Map.Entry<String,Integer> e : diCounts.entrySet()) {
            double freq = diWindows > 0 ? (double)e.getValue() / diWindows : 0.0;
            System.out.printf("%s\t%d\t%.6f%n", e.getKey(), e.getValue(), freq);
        }

        System.out.println();
        System.out.println("Observed trinucleotides (count, frequency = count/" + triWindows + "):");
        for (Map.Entry<String,Integer> e : triCounts.entrySet()) {
            double freq = triWindows > 0 ? (double)e.getValue() / triWindows : 0.0;
            System.out.printf("%s\t%d\t%.6f%n", e.getKey(), e.getValue(), freq);
        }
    }
}
