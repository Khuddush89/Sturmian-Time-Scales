# 🎨 Visual assets

| File / directory | Purpose |
| :--- | :--- |
| `hero.svg` | Gradient hero banner, responsive `viewBox`, accessible title |
| `typing.svg` | Local animated typing heading, with a static text fallback |
| `divider.svg` | Purple → pink → cyan section divider |
| `footer-wave.svg` | Animated decorative footer wave |
| `logo.svg` | Compact repository mark |
| `social-preview.svg`, `social-preview.png` | Editable source and 1280 × 640 GitHub preview |
| `demo-placeholder.svg` | Clearly marked terminal-recording placeholder |
| `screenshots/` | Seven actual research figures |

The README uses readme-typing-svg for its animated heading. Replace its external
image URL with `assets/typing.svg` for a completely local heading. SVG animations
are decorative; all files have useful static content.

## 📸 Additional screenshot slots

Add your own terminal capture as `screenshots/verification-terminal.png`
or a notebook capture as `screenshots/notebook-example.png`, then link it from
the README. These slots are placeholders, not claims that a notebook UI exists.

## 🎬 GIF slot

Record a real local run and add `demos/verification.gif`. Avoid embedding secrets
or personal terminal paths. Replace the demo placeholder with:

```markdown
![Verification and figure generation](assets/demos/verification.gif)
```

The color palette is purple `#8B5CF6`, pink `#EC4899`, and cyan `#06B6D4`.
Research plots retain their scientific colors and domain conventions.
