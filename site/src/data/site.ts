import {
          Activity, BookOpen, BrainCircuit, ClipboardList, Database, FileDown, Gauge, GitCompareArrows, Library, Milestone, Network, Radar, Route, Sparkles, Waves
        } from "lucide-react";

        const icons = {
          Activity, BookOpen, BrainCircuit, ClipboardList, Database, FileDown, Gauge, GitCompareArrows, Library, Milestone, Network, Radar, Route, Sparkles, Waves
        };

        export const paper = {
          title: "The Limits of Falsifiability",
          subtitle: "Dimensionality, measurement thresholds, and the sub-Landauer domain in biological systems",
          journal: "BioSystems",
          manuscriptId: "doi:10.1016/j.biosystems.2025.105608",
          decision: "Published",
          dueDate: "2025",
          doi: "10.1016/j.biosystems.2025.105608",
          abstract: "Falsifiability remains useful in biology only where the relevant signal is readable, the projection loss is limited, and the framework defining the test is explicit.",
          shortThesis: "Falsifiability remains useful in biology only where the relevant signal is readable, the projection loss is limited, and the framework defining the test is explicit.",
          mark: "F",
          pdfLabel: "PDF",
        };

        export const navLinks = [
          { href: "/", label: "Manuscript", icon: BookOpen },
          { href: "/dashboard", label: "Overview", icon: Radar },
          { href: "/revision", label: "Map", icon: ClipboardList },
          { href: "/audit", label: "Audit", icon: Gauge },
          { href: "/references", label: "References", icon: Library },
        ];

        export const downloads = [
  {
    "href": "/paper.pdf",
    "label": "PDF",
    "icon": "FileDown"
  },
  {
    "href": "/paper-source.tex",
    "label": "Source TeX",
    "icon": "FileDown"
  },
  {
    "href": "/citation-use-audit.md",
    "label": "Citation audit",
    "icon": "Database"
  },
  {
    "href": "/references.bib",
    "label": "Bibliography",
    "icon": "Library"
  }
].map((item) => ({ ...item, icon: icons[item.icon as keyof typeof icons] }));
        export const headlineMetrics = [
  {
    "label": "Publication",
    "value": "BioSystems 258",
    "detail": "Article 105608 (2025)",
    "tone": "primary"
  },
  {
    "label": "Paper DOI",
    "value": "105608",
    "detail": "10.1016/j.biosystems.2025.105608",
    "tone": "tertiary"
  },
  {
    "label": "Cited sources",
    "value": "16",
    "detail": "16 bibliography entries parsed from TeX",
    "tone": "secondary"
  },
  {
    "label": "PDF/text coverage",
    "value": "14/16",
    "detail": "0 DOI-backed entries still need PaperLibrary harvest",
    "tone": "success"
  }
];
        export const stanceCards = [
  {
    "title": "Framework as projection",
    "body": "Falsifiability depends on the measurement, representation, and decision maps that define the framework before evidence is even gathered.",
    "icon": "Milestone"
  },
  {
    "title": "Physical readability",
    "body": "Sub-Landauer biological signals can be causally potent without being reliably available as single-event binary records.",
    "icon": "Gauge"
  },
  {
    "title": "Projection loss",
    "body": "High-dimensional biological states can survive aggregate measurement while resisting binary single-trial falsification.",
    "icon": "Network"
  },
  {
    "title": "Ensemble epistemology",
    "body": "The stronger biological move is not abandoning falsification, but moving scale-specific claims into ensemble and multi-scale inference.",
    "icon": "Radar"
  }
].map((item) => ({ ...item, icon: icons[item.icon as keyof typeof icons] }));
        export const revisionPriorities = [
  {
    "title": "Separate three limitations",
    "reviewer": "Paper architecture",
    "priority": "Critical",
    "body": "Physical measurement, dimensional projection, and framework dependence are kept analytically distinct."
  },
  {
    "title": "Use Wigner as selection bias",
    "reviewer": "Paper architecture",
    "priority": "High",
    "body": "The effectiveness of mathematics is recast as domain selection for low projection-loss systems."
  },
  {
    "title": "Preserve falsification",
    "reviewer": "Paper architecture",
    "priority": "High",
    "body": "The paper does not reject Popper; it localizes where binary epistemology is strongest."
  },
  {
    "title": "Connect to timing limits",
    "reviewer": "Paper architecture",
    "priority": "Medium",
    "body": "The companion timing paper is used as the mechanistic extension of framework-limited falsifiability."
  }
];
        export const stressTests = [
  {
    "title": "Framework as projection",
    "body": "Falsifiability depends on the measurement, representation, and decision maps that define the framework before evidence is even gathered.",
    "icon": "Milestone"
  },
  {
    "title": "Physical readability",
    "body": "Sub-Landauer biological signals can be causally potent without being reliably available as single-event binary records.",
    "icon": "Gauge"
  },
  {
    "title": "Projection loss",
    "body": "High-dimensional biological states can survive aggregate measurement while resisting binary single-trial falsification.",
    "icon": "Network"
  },
  {
    "title": "Ensemble epistemology",
    "body": "The stronger biological move is not abandoning falsification, but moving scale-specific claims into ensemble and multi-scale inference.",
    "icon": "Radar"
  }
].map((item) => ({ ...item, icon: icons[item.icon as keyof typeof icons] }));
        export const empiricalTests = [
  {
    "title": "Readable signal",
    "body": "Can the relevant biological state be stabilized as a reliable record at the scale being tested?"
  },
  {
    "title": "Projection budget",
    "body": "Does the measurement preserve the causal distinctions the hypothesis actually depends on?"
  },
  {
    "title": "Framework agreement",
    "body": "Are observers using the same measurement, representation, and decision boundary?"
  }
];
        export const evidenceClusters = [
  {
    "title": "Epistemic foundations",
    "keys": [
      "popper1959",
      "duhem1906",
      "quine1951"
    ],
    "icon": "Library"
  },
  {
    "title": "Physical measurement limits",
    "keys": [
      "landauer1961",
      "engel2007",
      "cao2020"
    ],
    "icon": "Gauge"
  },
  {
    "title": "Sub-Landauer biological effects",
    "keys": [
      "anastassiou2011",
      "chiang2019",
      "mcdonnell2011",
      "stocks2000"
    ],
    "icon": "Waves"
  },
  {
    "title": "Dimensional projection",
    "keys": [
      "wigner1960",
      "dill2012",
      "todd2025timing"
    ],
    "icon": "Network"
  }
].map((item) => ({ ...item, icon: icons[item.icon as keyof typeof icons] }));
        export const reviewerPosture = [
  {
    "reviewer": "Published article",
    "stance": "BioSystems",
    "summary": "This native site renders the full published manuscript for Limits of Falsifiability as web text with live citations and local PDF exports."
  },
  {
    "reviewer": "Reference layer",
    "stance": "PaperLibrary-backed",
    "summary": "Each cited key receives a reference page with manuscript contexts, DOI links where available, and the current local PDF/text harvest state."
  }
];
        export const reframes = [
  {
    "before": "A hypothesis is falsifiable or not.",
    "after": "A hypothesis is falsifiable relative to the framework that projects states into evidence."
  },
  {
    "before": "Measurement only observes biology.",
    "after": "Measurement also selects and compresses the biological state space."
  },
  {
    "before": "Mathematical success is universal.",
    "after": "Mathematical success is strongest where projection loss is unusually low."
  }
];
        export const detailHighlights = [
  {
    "id": "main-thesis",
    "phrase": "Falsifiability remains useful in biology only where the relevant signal is readable, the projection loss is limited, and the framework defining the test is explicit.",
    "title": "Main thesis",
    "summary": "Falsifiability remains useful in biology only where the relevant signal is readable, the projection loss is limited, and the framework defining the test is explicit.",
    "bullets": [
      "Falsifiability depends on the measurement, representation, and decision maps that define the framework before evidence is even gathered.",
      "Sub-Landauer biological signals can be causally potent without being reliably available as single-event binary records.",
      "High-dimensional biological states can survive aggregate measurement while resisting binary single-trial falsification."
    ]
  }
];
