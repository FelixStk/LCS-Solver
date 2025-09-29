int calculateTotal(const std::span<int> items) {
  int sum = 0;
  for (int number : items) {
    sum += number;
  }
  return sum;
}