# Claude Code — Project Notes
Do not write any code until I give the go ahead. Give me an estimate of the amount of code you intend to write. 

## Working directory
Always work directly in this repository on the `main` branch.
Do not create worktrees or feature branches unless explicitly asked.

## Color rules
All colors must come from `assets/css/palette.css` CSS custom properties.
Never write a raw hex value (#rrggbb) in any applet file.
See `Writing_Applet.txt` → COLORS section for the full guide.

## Applet files
All applet HTML partials live in `_includes/applets/`.
Jekyll includes use:  {% include applets/mymodel_applet.html %}
