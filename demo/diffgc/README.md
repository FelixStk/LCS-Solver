# Demo: diffgc
Run a lcs algorithm to generate the difference between two files.
```shell 
diffgc res/file1.cc res/file2.cc -i
```

## Usage
```text
Usage: diffgc.exe file1 file2 [-j config.json] [-s] [--gnu]
Parameters:
file1 file2    : Filepaths to the original and modified file
Options:
-j config.json : A JSON file that contains zero or one constraint
-s file.json   : Save a constraint to a JSON file
-i             : Prompt for a constraint type and parameters
--gnu          : Print in the unified output format
--help or -h   : Display this usage message
--version      : Display the version of the program
```

## Notes:
- The `--gnu` flag specifies that the difference should be displayed in [the unified diff format](https://www.gnu.org/software/diffutils/manual/html_node/Unified-Format.html). If this flag is not set, the output format defaults to a standard one-window view of the differences.
- If the save flag `-s` is set (in combination with `-i`), the first (and only) constraint is saved as JSON, keyed with `constraint`, using `lcs_solver::constraints::to_json`.
- Every line in `file1` and `file2` is mapped to a symbol by `diffgc::FileAdapter.cc`. The string type used in `lcs_solver` is hardcoded in `include/lcs_solver/util/CommonTypes.h`. This demo uses `Symbol = char32_t`.

## Example:
Assume:
- `res/file1.cc`: 
    ```c++
    int calculateTotal(const std::span<int> items) {
      int sum = 0;
      for (int number : items) {
        sum += number;
      }
      return sum;
    }
    ```
- `res/file2.cc`
    ```c++
    #include <span>
    
    int calculateTotal(const std::span<int> items) {
      int sum = 0;
      if (items.empty()) return sum;
      sum = items.front() + calculateTotal(items.subspan(1));
      return sum;
    }
    ```
- A relaxed MC constraint in `config_mc.json`: 
    ```json
    {
      "constraint": {
        "data": {
          "gap_": [
            [0, 8],
            [0, 8],
            [0, 8],
            [0, 8],
            [0, 8],
            [0, 8]
          ],
          "gap_lower_bound_": 0,
          "gap_max_length_": 8,
          "gap_num_uniques_": 1,
          "gap_upper_bound_": 8,
          "t_": 6
        },
        "identifier": "Constraint_MC"
      },
      "adapter": {
        "Files": [
          "res/file1.cc",
          "res/file2.cc"
        ],
        "MapR": [
          [107, "  sum = items.front() + calculateTotal(items.subspan(1));"],
          [106, "  if (items.empty()) return sum;"],
          [105, ""],
          [104, "#include <span>"],
          [103, "}"],
          [102, "  return sum;"],
          [101, "  }"],
          [100, "    sum += number;"],
          [99, "  for (int number : items) {"],
          [98, "  int sum = 0;"],
          [97, "int calculateTotal(const std::span<int> items) {"]
        ]
      }
    }
    ```
Then the output of `diffgc res/file1.cc res/file2.cc -j res/config_mc.json` will be:
```text
+ | #include <span>
+ |
= | int calculateTotal(const std::span<int> items) {
= |   int sum = 0;
- |   for (int number : items) {
- |     sum += number;
- |   }
+ |   if (items.empty()) return sum;
+ |   sum = items.front() + calculateTotal(items.subspan(1));
= |   return sum;
= | }
```
And the output of `diffgc res/file1.cc res/file2.cc -j res/config_mc.json --gnu` will something like this:
```text
--- res/file1.cc        YYYY-MM-DD HH:MM:SS.nnnnnnnnn +/-HHMM
+++ res/file2.cc        YYYY-MM-DD HH:MM:SS.nnnnnnnnn +/-HHMM
@@ -1,0 +1,2 @@
+#include <span>
+
@@ -3,3 +5,2 @@
-  for (int number : items) {
-    sum += number;
-  }
+  if (items.empty()) return sum;
+  sum = items.front() + calculateTotal(items.subspan(1));
```