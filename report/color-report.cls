\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{color-report}[2021/06/01 colorful report class]

% Load the base document class
\LoadClassWithOptions{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PACKAGE LOADING                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphics
\RequirePackage{geometry}		% For page layout
\RequirePackage{tikz}           % Nice drawings
\usetikzlibrary{
    shapes.geometric, 			% For geometric shapes (rectangle ...)
	calc, 						% For coordinates operations
	positioning					% For relative positioning
    }
\RequirePackage{hyperref}		% For better links in document


% Maths
\RequirePackage{amsmath, amssymb}   % For math symbols and commands

% Special math environments
% NOTE : must be included after babel, amsmath and hyperref and should be used without
% asmthm
\RequirePackage[amsmath, hyperref]{ntheorem}	% For nice theorem-like environnements (replaces amsthm)
\RequirePackage[ntheorem]{mdframed}				% Create framed boxes
\RequirePackage[ruled,vlined]{algorithm2e}		% For pseudocode
\RequirePackage[outputdir=build]{minted}		% For code blocks

% Latex programming
\RequirePackage{etoolbox}			% For control flow helpers

% Fonts
\RequirePackage{fontspec}			% Easily load fonts
\RequirePackage{helvet}				% For the helvetica font (similar to Arial)
% \RequirePackage[
% 	rm, 
% 	mono=false, 
% 	type1=true]{libertine}			% Use the libertine font as the roman font
% \RequirePackage{unicode-math}		% For easily changing the maths fonts	
% \setmathfont[
% 	Scale=MatchUppercase
% ]{libertinusmath-regular.otf}		% Set the math font to libertinus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title and information
\NewDocumentCommand{\subtitle}{m}{\NewDocumentCommand{\@subtitle}{}{#1}}
\NewDocumentCommand{\reportType}{m}{\NewDocumentCommand{\@reportType}{}{#1}}
\NewDocumentCommand{\@info}{}{}
\NewDocumentCommand{\info}{m m}{\g@addto@macro\@info{%
	\vspace{1em}
	{\Large\bfseries \textit{#1} \par}
	{%
		\large%
		#2 \par
	}
}}

% Logos
\NewDocumentCommand{\mainLogo}{m}{\NewDocumentCommand{\@mainLogo}{}{#1}}
\NewDocumentCommand{\secondaryLogo}{m}{\NewDocumentCommand{\@secondaryLogo}{}{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STYLES                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default colors
\definecolor{black}{RGB}{0, 0, 0}
\definecolor{black-light}{RGB}{44, 46, 53}
\definecolor{primary}{RGB}{0, 150, 170}
\definecolor{accent-1}{RGB}{255, 125, 40}
\definecolor{accent-2}{RGB}{0, 160, 130}
\definecolor{accent-2-light}{RGB}{140, 190, 170}
\definecolor{accent-3}{RGB}{250, 100, 230}
\definecolor{accent-3-light}{RGB}{200, 160, 190}

% itemize color and marker
\RenewDocumentCommand{\labelitemi}{}{\color{accent-1}$\bullet$}
\RenewDocumentCommand{\labelitemii}{}{\color{accent-1}$\ast$}
\RenewDocumentCommand{\labelitemiii}{}{\color{accent-1}$\diamond$}
\RenewDocumentCommand{\labelitemiv}{}{\color{accent-1}$\cdot$}

% links
\hypersetup{
	colorlinks=true,
	linkcolor=black-light,
	citecolor=accent-1,
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENVIRONMENTS                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\theoremstyle{plain}

% Remark
\theoremheaderfont{\color{accent-2}\bfseries}
\newmdtheoremenv[
	linewidth=5pt, 
	linecolor=accent-2-light, 
	topline=false, 
	bottomline=false, 
	rightline=false,
]{remark}{Remark}

% Important proposition
\theoremheaderfont{\color{primary}\bfseries}
\newmdtheoremenv[
	linewidth=5pt, 
	linecolor=primary-light, 
	topline=false, 
	bottomline=false, 
	rightline=false,
]{iproposition}{Proposition}

\theoremheaderfont{\color{black}\bfseries}
\newtheorem{proposition}{Proposition}
\newtheorem{definition}{Definition}

\theoremstyle{nonumberplain}
\newtheorem{proof}{Proof}

% algorithm2e environment
\SetKw{Continue}{continue}
\SetKwInput{KwStaticVariables}{Static variables}

% minted environment
\BeforeBeginEnvironment{minted}{%
	\begin{mdframed}[linewidth=5pt, 
		linecolor=accent-3-light, 
		topline=false, 
		bottomline=false, 
		rightline=false]
}
\AfterEndEnvironment{minted}{\end{mdframed}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE PAGE                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RenewDocumentCommand {\maketitle} {} {%
	\newgeometry{left=0pt, right=0pt, top=0pt, bottom=0pt}
	\begin{tikzpicture}[remember picture, overlay]
		% The blue vertical banner
		\fill[primary] (current page.north west) rectangle ($(current page.south west) + (25em, 0)$);

		% The report type
		\node[anchor=north west, rotate=90] at ($(current page.south west) + (2em, 2em)$) {
			{\fontsize{41}{41}\selectfont\color{white}\sffamily\bfseries \@reportType}
		};
		
		% The main logo
		\node[anchor=north west, inner sep=0, outer sep=0] at ($(current page.north west) + (0, 0)$) {
			\includegraphics[width=23em]{\@mainLogo}
		};

		% The secondary logo
		\node[anchor=north west, inner sep=0, outer sep=0] at ($(current page.north west) + (25em + 2.5em, 0)$) {
			\includegraphics[width=23em]{\@secondaryLogo}
		};

		% The title
		\node[anchor=north west, minimum height=3cm] (title) at ($(current page.north west) + (25em + 2em, -20em)$) {
			\begin{minipage}[r]{.5\textwidth}
				{\flushright\color{black-light}\Huge\selectfont\sffamily\bfseries \@title \par}
				{\flushright\LARGE\selectfont \textit{\@subtitle} \par}
			\end{minipage}
		};

		% The author(s)
		\node[anchor=north, below=of title, inner sep=0, outer sep=2em] (author) {
			\begin{minipage}{.5\textwidth}
				{\centering\huge\selectfont\sffamily\bfseries \@author \par}
			\end{minipage}
		};

		% The date
		\node[anchor=north, below=of author, inner sep=0, outer sep=0] {
			{\Large\selectfont \@date \par}
		};
		
		% The information fields
		\node[anchor=south west, minimum height=3cm] at ($(current page.south west) + (25em + 2em, 2em)$) {		
			\begin{minipage}{.5\textwidth}
				\@info
			\end{minipage}
		};
	\end{tikzpicture}
    \restoregeometry{}
}