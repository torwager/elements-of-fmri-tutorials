# Elements of fMRI Analysis — Interactive Tutorials

Hands-on tutorials accompanying **_Elements of Functional Magnetic Resonance Imaging_** by Tor D. Wager and Martin A. Lindquist.

- 🌐 **Site:** https://torwager.github.io/elements-of-fmri-tutorials/ *(goes live at launch)*
- 📖 42 chapters across 7 parts, each with concept summaries, runnable MATLAB + Python code, thought questions, and self-quizzes
- 🧰 MATLAB exercises use [CANlab tools](https://canlab.github.io); Python exercises use [nilearn](https://nilearn.github.io) and friends
- ▶️ Much of the Python code runs in-browser (WebAssembly) — no installation needed

## Development

Built with [Jupyter Book 2 / MyST](https://mystmd.org).

```bash
npm install -g mystmd
myst start          # local preview of the book
```

The custom landing page lives in `landing/` and is assembled with the built book by the deploy workflow (`.github/workflows/deploy.yml`): landing page at the site root, book under `/book/`.

## License

Content: CC-BY-4.0 · Code: MIT
