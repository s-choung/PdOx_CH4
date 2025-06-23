\section*{PdOx\_CH\textsubscript{4}}

\subsection*{1. Scripts}

All automation tools for generating input and analyzing NEB results are provided in the \texttt{scripts/} folder:

\begin{itemize}
    \item \texttt{generate\_initial\_final\_states.py}: Places CH\textsubscript{4} molecule and generates CH\textsubscript{3} + H final geometries with overlap checks.
    \item \texttt{interpolate\_neb\_images.py}: Generates NEB images between valid initial/final state pairs.
    \item \texttt{run\_neb\_job.sh}: Shell script template to submit NEB jobs to a cluster or local environment.
    \item \texttt{visualize\_trajectories.py}: Generates overlayed CH\textsubscript{4} activation geometries for each model.
\end{itemize}

\subsection*{2. Optimized Surface Structures}

The \texttt{input\_structures/} folder contains optimized surface models (in \texttt{.vasp} format) used to generate activation pathways:

\begin{itemize}
    \item \texttt{Pd.vasp}, \texttt{PdO.vasp}, \texttt{Pd-PdO.vasp}, and \texttt{PdOx.vasp}
\end{itemize}

\subsection*{3. NEB Output Structures}

The \texttt{neb\_outputs/} directory contains:

\begin{itemize}
    \item Initial and final geometries for each surface site.
    \item Final NEB-converged transition states.
    \item Site-wise organization by surface type and site index.
\end{itemize}

\subsection*{4. Visualization}

Overlay plots of all CH\textsubscript{4} activation trajectories are available in the \texttt{visualizations/} folder:

\begin{itemize}
    \item These \texttt{.xyz} files can be opened in visualization tools (e.g., VESTA, ASE GUI, OVITO) to examine the distribution of reaction geometries across all surface sites.
    \item Color code: Pd (turquoise), O (red), Ce (beige), C (gray), H (white).
\end{itemize}

\subsection*{Citation}

If you use this repository or scripts in your work, please cite:

\begin{quote}
    \textbf{[DOI Placeholder]} ``Title of Manuscript'', \textit{Journal}, 2025.
\end{quote}

\noindent For any questions or issues, feel free to open an issue or contact the authors.
