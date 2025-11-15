Make this diff AS ULTRA-MINIMAL as possible (compared to main) while retaining functional equivalence.

Your job is to look at the git diff of this branch compared to main, and remove
AS MUCH CODE AS POSSIBLE while ensuring that the resulting code still behaves
the same way as before.  The goal is to allow a human reader to quickly
understand the core change.

DO:

- remove all but the most important inline comments
- leave only a one-line docstring at the top of each function or class, with
  only the most essential information
- eliminate verbosity, duplication, and unnecessary code
- simplify logic where possible, using more concise constructs,
  list comprehensions, assignment expressions, built-in functions, etc.
- remove overly defensive code that is not strictly necessary for correct operation

DO NOT:

- change functionality.
- add new functionality.
- obfuscate the code or use clever tricks that reduce readability.
- remove code that existed on main.

IMPORTANT:

- Any tests that pass before you begin should still pass after you finish.
- You should NEVER run `git commit` or `git push` commands yourself.
  Just modify the code to be ultra-minimal and stop there.
