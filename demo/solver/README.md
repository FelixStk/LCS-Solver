# Demo: solver
This demo showcases the library's functionality for generating, processing, and solving LCS problems

## Quickstart
The program has two modes: If the user calls the program with the `-i` flag, it will be started in the _interactive mode_. It allows the user to interactively create, load, save, print and solve lcs problems. The other mode of the program is the _file mode_. If the user calls the program with `-f problem.json`, the solver will try to parse the problem and solve it.

For usage examples, see the integration tests in `test_solver.cc`. Running a specific test will generate a problem-type-specific `test.json` file in the build directory (where the LCS_Test executable is located).

## Interactive mode
In the interactive mode, the user controls the actions of the program by entering characters in a main loop:
- `n` to create a new problem (calls `ProblemFactory::CreateFromDialog()`)
- `o` to open a problem from a JSON file (calls `ProblemFactory::CreateFromJson(file_name)`)
- `s` to save a problem to a JSON file (calls `ProblemFactory::SaveToJson(file_name, problem_ptr)`)
- `x` to execute an algorithm defined in the problem generation process (calls `problem_ptr->ExecuteAlgo()`)
- `p` to print the currently loaded or specified problem to the standard output
- `e` to exit the action loop

During the creation of a problem, the user must choose a type of problem to generate. `ProblemType::LCS_Base` can be used to specify a general problem by composing arbitrary many strings and constraints. If it is chosen, there is no safeguard in place that ensures that the problem is well-defined. The other `ProblemTypes` are introduced by Adamson et al. in [1]. 

## File mode
In the file mode, a JSON is parsed to define an LCS Problem. The JSON file must include valid values for the type, name, description, strings, constraint map and algorithm. Here is an example of a valid JSON:
```json
{
  "type_": "LCS_Sigma",
  "name_": "my_problem_name",
  "description_": "my_problem_desc",
  "algo_": "LLCS2_SA_RMQ",
  "spv_": [
    "test",
    "test"
  ],
  "map_": [
    [
      {"ConstraintType": "LCS_2"},
      {
        "data": {"k_": 2},
        "identifier": "LCS_2"
      }
    ],
    [
      {"ConstraintType": "LCS_Const_Sigma"},
      {
        "data": {"sigma_": 3},
        "identifier": "LCS_Const_Sigma"
      }
    ],
    [
      {"ConstraintType": "Constraint_Sigma"},
      {
        "data": {
          "left_": [
            [116, [0, 4]],
            [115, [0, 4]],
            [101, [0, 4]]],
          "right_": [
            [116, [0, 4]],
            [115, [0, 4]],
            [101, [0, 4]]
          ]
        },
        "identifier": "Constraint_Sigma"
      }
    ]
  ]
}
```
For more information, please refer to the respective header files of the objects:
- `type_` is string is defined by the Macro in `include/lcs_solver/problems/ProblemType.h`
- `name_` and `description_` are normal utf8 strings
- `spv_` is used for initializing a vector of string pointer. In the JSON file it is an array of strings.
- `algo_` is string is defined by the Macro in `include/lcs_solver/algorithms/AlgoType.h`
- `map_` is an array of key-value pairs to initialize a `ConstraintMap`. The string to identify a constraint is defined in the field `kName` of a constraint. The value can be accessed by the function using the map getter `GetMapTypeToName` declared in `ConstraintType.h`.

## References

[1] Adamson, Duncan et al. "Longest Common Subsequence with Gap Constraints" arXiv, 2023.