import java.util.*;
import java.util.stream.*;

public class Lab1_1 {
    public static String alphabetOf(String s) {
        return String.join("", s.chars()
                .mapToObj(c -> (char) c)
                .map(Object::toString)
                .collect(Collectors.toCollection((TreeSet::new))));
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.print("Enter input: ");
        String s = scanner.nextLine();
        System.out.println(alphabetOf(s));
        scanner.close();
    }
}
