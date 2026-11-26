% Notation macros: HTML gets sub/superscripts, text help and the PDF manual
% get the plain ASCII spelling (which is also the output column name).
\newcommand{\sb}{\ifelse{html}{\out{#1<sub>#2</sub>}}{#1_#2}}
\newcommand{\sps}{\ifelse{html}{\out{#1<sup>#2</sup>}}{#1^#2}}
\newcommand{\inv}{\ifelse{html}{\out{#1<sup>&minus;1</sup>}}{#1^-1}}
\newcommand{\hsq}{\ifelse{html}{\out{h&sup2;}}{h2}}
\newcommand{\sqrthsq}{\ifelse{html}{\out{&radic;h&sup2;}}{sqrt(h2)}}
\newcommand{\GxE}{\ifelse{html}{\out{G&times;E}}{GxE}}
\newcommand{\VGx}{\ifelse{html}{\out{V<sub>G&times;#1</sub>}}{V_Gx#1}}
\newcommand{\by}{\ifelse{html}{\out{&times;}}{x}}
\newcommand{\rt}{\ifelse{html}{\out{&radic;(}#1\out{)}}{sqrt(#1)}}
