# -----------------------------
# Options affecting formatting.
# -----------------------------
with section("format"):
    # How wide to allow formatted cmake files
    line_width = 100

    # How many spaces to tab for indent
    tab_size = 4

    # If an argument group contains more than this many sub-groups (parg or kwarg
    # groups) then force it to a vertical layout.
    max_subgroups_hwrap = 2

    # If a positional argument group contains more than this many arguments, then
    # force it to a vertical layout.
    max_pargs_hwrap = 4

    # If a cmdline positional group consumes more than this many lines without
    # nesting, then invalidate the layout (and nest)
    max_rows_cmdline = 2

    # If a statement is wrapped to more than one line, than dangle the closing
    # parenthesis on its own line.
    dangle_parens = True

    # If the trailing parenthesis must be 'dangled' on its on line, then align it
    # to this reference: `prefix`: the start of the statement,  `prefix-indent`:
    # the start of the statement, plus one indentation  level, `child`: align to
    # the column of the arguments
    dangle_align = "prefix"

    # If a candidate layout is wrapped horizontally but it exceeds this many
    # lines, then reject the layout.
    max_lines_hwrap = 2

    # What style line endings to use in the output.
    line_ending = "unix"

    # Format command names consistently as 'lower' or 'upper' case
    command_case = "lower"

    # Format keywords consistently as 'lower' or 'upper' case
    keyword_case = "upper"

    # If true, the argument lists which are known to be sortable will be sorted
    # lexicographicall
    enable_sort = True

    # If true, the parsers may infer whether or not an argument list is sortable
    # (without annotation).
    autosort = True

    # By default, if cmake-format cannot successfully fit everything into the
    # desired linewidth it will apply the last, most aggressive attempt that it
    # made. If this flag is True, however, cmake-format will print error, exit
    # with non-zero status code, and write-out nothing
    require_valid_layout = False

# ----------------------------
# Options affecting the linter
# ----------------------------
with section("lint"):
    # a list of lint codes to disable
    disabled_codes = ["C0103"]
