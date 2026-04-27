# Claude Code — Project Notes

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
