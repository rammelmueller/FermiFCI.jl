The script `single_diag.jl` . It can be invoked as follows (from within this directory):

```
    julia single_diag.jl
```
or, after starting the REPL with `julia`, the script may simply be included via
```
    julia> include("single_diag.jl")
```
Note that the second way circumvents the need of JIT compilation every time the script is invoked and therefore performs *much* faster when called the second time.

Input parameters may be adjusted in the YAML file `config.yml`.
