import numpy as np


class DiscreteStepPredictor:
    def __init__(self, matrix, initial_vector):
        self.matrix = np.array(matrix)
        self.state = np.array(initial_vector)
        self.initial_state = np.array(initial_vector)

        self._validate_inputs()

    def _validate_inputs(self):
        if len(self.matrix.shape) != 2:
            raise ValueError("Matrix must be 2-dimensional")

        if self.matrix.shape[0] != self.matrix.shape[1]:
            raise ValueError("Matrix must be square")

        if self.matrix.shape[1] != len(self.state):
            raise ValueError("Matrix columns must match vector length")

    def reset(self):
        self.state = self.initial_state.copy()

    def predict_steps(self, num_steps=5):
        predictions = [self.state.copy()]

        for step in range(1, num_steps + 1):
            self.state = self.matrix @ self.state
            predictions.append(self.state.copy())

        return predictions

    def get_state_probabilities(self):
        if np.all(self.matrix >= 0) and np.allclose(self.matrix.sum(axis=1), 1):
            return self.state / self.state.sum()
        return self.state

    def print_predictions(self, num_steps=5, precision=4):
        predictions = self.predict_steps(num_steps)

        print(f"DISCRETE STEP PREDICTION - {num_steps} STEPS")
        print(f"Matrix size: {self.matrix.shape[0]}x{self.matrix.shape[1]}")
        print(f"Initial vector: {self.initial_state}")

        for i, state in enumerate(predictions):
            if i == 0:
                print(f"Step {i} (Initial): {np.round(state, precision)}")
            else:
                print(f"Step {i}:         {np.round(state, precision)}")


def main():
    while True:
        print("\n" + "=" * 60)
        print("DISCRETE STEP PREDICTION SYSTEM")
        print("=" * 60)
        print("1. Enter custom matrix and vector")
        print("2. Exit")

        choice = input("Select option: ").strip()

        if choice == '1':
            try:
                n = int(input("\nEnter matrix size (n for nxn): "))

                print(f"\nEnter {n}x{n} matrix:")
                matrix = []
                for i in range(n):
                    while True:
                        try:
                            row_input = input(f"Row {i + 1} (enter {n} space-separated numbers): ")
                            row = list(map(float, row_input.split()))
                            if len(row) != n:
                                print(f"Error: Expected {n} numbers")
                                continue
                            matrix.append(row)
                            break
                        except ValueError:
                            print("Error: Please enter valid numbers")

                print(f"\nEnter initial vector ({n} values):")
                while True:
                    try:
                        vector_input = input("Enter space-separated values: ")
                        vector = list(map(float, vector_input.split()))
                        if len(vector) != n:
                            print(f"Error: Expected {n} values")
                            continue
                        break
                    except ValueError:
                        print("Error: Please enter valid numbers")

                predictor = DiscreteStepPredictor(matrix, vector)

                while True:
                    print("PREDICTION OPTIONS")
                    steps = int(input("Enter number of steps to predict: "))
                    predictor.print_predictions(steps)

                    cont = input("\nMake another prediction with same system? (y/n): ").lower()
                    if cont != 'y':
                        predictor.reset()
                        break

            except Exception as e:
                print(f"Error: {e}")

        elif choice == '2':
            print("Exiting program.")
            break
        else:
            print("Invalid choice. Please try again.")


if __name__ == "__main__":
    main()
