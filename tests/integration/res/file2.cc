#include <span>

int calculateTotal(const std::span<int> items) {
  int sum = 0;
  if (items.empty()) return sum;
  sum = items.front() + calculateTotal(items.subspan(1));
  return sum;
}