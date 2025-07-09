
Take inventory of boltz code base
make no to limit changes to source code -- build funtcions that integrate into/on top of code base
ASK when chaning boltz code
create / look at  baseline.md file to remind you what you observe


Summarize progress in a PROGRESS.md file
check plan.md for implementation

NO interface{} or any{} - use concrete types!
NO time.Sleep() or busy waits - use channels for synchronization!
NO keeping old and new code together
NO migration functions or compatibility layers
NO versioned function names (processV2, handleNew)
NO custom error struct hierarchies
NO TODOs in final code
Meaningful names: userID not id

Stop - Don't spiral into complex solutions
Step back - Re-read the requirements --> look at plan.md
Simplify - The simple solution is usually correct
Ask - "I see two approaches: [A] vs [B]. Which do you prefer?"


DELETE old code