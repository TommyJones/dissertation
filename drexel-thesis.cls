%% This file was generated from the drexel-thesis
%% LaTeX class and then modified for use with the
%% dragondown package. If this file is to be replaced
%% at any point with an updated version of the
%% drexel-thesis class, please note the differences
%% between this file and the one generated to ensure
%% that no breaking changed will occur.
%%
%%
%% This is file `drexel-thesis.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% drexel-thesis.dtx  (with options: `class')
%%
%% This is a generated file.
%%
%% Copyright (C) 2010-2012 W. Trevor King
%%
%% This file may be distributed and/or modified under the conditions of
%% the LaTeX Project Public License, either version 1.3 of this license
%% or (at your option) any later version.  The latest version of this
%% license is in:
%%
%%    http://www.latex-project.org/lppl.txt
%%
%% and version 1.3 or later is part of all distributions of LaTeX version
%% 2003/06/01 or later.
%%
\def\fileversion{0.15}
\def\filedate{2013/05/21}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{drexel-thesis}
 [\filedate\ \fileversion\ Drexel Thesis]
\RequirePackage{xkeyval}
\newif\if@mainmatter \@mainmattertrue
\newif\if@DUT@subfigopt
\newif\if@DUT@pagebreak
\newif\if@dottedchapter
\newif\if@boldchapter
\newif\if@boldpart
\newif\if@dottedpart
\newif\if@indenttoc
\newif\if@blackref
\newif\if@draftwatermark
\newif\if@xdraftwatermark
\newif\if@indentfirst
\newif\if@final
\newif\if@daring
\newif\if@approvalform
\newcommand\iffinal[2]{\if@final #1\else #2\fi}
\let\DUTiffinal\iffinal
\newcommand\ifdraftwatermark[2]{\if@draftwatermark #1\else #2\fi}
\let\DUTifdraftwatermark\ifdraftwatermark
\newcommand\ifdaring[2]{\if@daring #1\else #2\fi}
\let\DUTifdaring\ifdaring
\newcommand\DUT@makeapprovalform{
  % \includepdf[pages=-]{approval-form.pdf}
  \includepdf[pages=-]{\@approvalForm}
}
% \newcommand\DUT@makeapprovalform{
%   \thispagestyle{empty}
%   \begingroup
%   \newlength{\tempskipa}
%   \newlength{\tempskipb}
%   \parindent \z@
%   \settowidth{\tempskipa}{\drexellogobw}
%   \addtolength{\tempskipa}{1em}
%   \setlength{\tempskipb}{\textwidth}
%   \addtolength{\tempskipb}{-\tempskipa}
%   \addtolength{\tempskipb}{-1em}
%   \begin{minipage}[c]{\tempskipa}
%     \raisebox{-0.5\height}{\drexellogobw}
%   \end{minipage}
%   \hfill
%   \begin{minipage}[c]{\tempskipb}
%     \textbf{\Huge Office of Graduate Studies} \linebreak
%     \vskip 0.1\baselineskip
%     \textbf{\LARGE Dissertation / Thesis Approval Form} \linebreak
%   \end{minipage}
%
%   \vskip 2\baselineskip
%
%   This form is for use by all doctoral and master's students with a
%   dissertation/thesis requirement.  Please print clearly as the
%   library will bind a copy of this form with each copy of the
%   dissertation/thesis.  All doctoral dissertations must conform to
%   university format requirements, which is the responsibility of the
%   student and supervising professor.  Students should obtain a copy of
%   the Thesis Manual located on the library website.
%
%   \vspace*{\fill}
%   \settowidth{\tempskipb}{Dissertation/Thesis Title:}
%   \setlength{\tempskipa}{\textwidth}
%   \addtolength{\tempskipa}{-\tempskipb}
%   \addtolength{\tempskipa}{-4em}
%   \begin{spacing}{2}
%     \begin{tabular*}{\textwidth}{@{\extracolsep{\fill}} l p{\tempskipa}}
%       \textbf{Dissertation/Thesis Title:} & \@DUT@title \\
%       \textbf{Author:} & \@DUT@author
%     \end{tabular*}
%
%     \vspace*{\fill}
%
%     \begin{center}
%       \textbf{This dissertation/thesis is hereby accepted and approved.}
%     \end{center}
%
%     \vspace*{\fill}
%
%     \textbf{Signatures:} \\
%     Examining Committee \\
%     \begin{tabular*}{\textwidth}{@{\extracolsep{\fill}} p{0pt} l p{\tempskipa}}
%       & Chair            & \hrulefill \\
%       & Members          & \hrulefill \\
%       &                  & \hrulefill \\
%       &                  & \hrulefill \\
%       &                  & \hrulefill \\
%       &                  & \hrulefill \\
%       & Academic Advisor & \hrulefill \\
%       & Department Head  & \hrulefill
%     \end{tabular*}
%
%     \vspace*{\fill}
%   \end{spacing}
%
%   \begin{center}\scriptsize
%     Office of Graduate Studies
%     $\bullet$ 3141 Chestnut St.
%     $\bullet$ Randell Hall 240
%     $\bullet$ Philadelphia, PA 19104 \\
%     \textbf{Tel.}: 215-895-0366
%     $\bullet$ \textbf{Fax}: 215-895-0495
%     $\bullet$ \textbf{Email}:
%       \href{mailto:graduate@drexel.edu}{graduate@drexel.edu}
%     $\bullet$ \textbf{Web}:
%       \href{http://www.drexel.edu/provost/graduatestudies/}{%
%                    www.drexel.edu/provost/graduatestudies}
%   \end{center}
%   \endgroup
%   \pagebreak
% }

\DeclareOptionX{subfig}{\@DUT@subfigopttrue}
\DeclareOptionX{finalspace}{
  \setlength\overfullrule{0pt}
  \def\@DUT@toc@sep{1} % content page separation (# blanklines between entries)
  \def\@ssp{1} % single space
  \def\@dsp{2} % double space
  \def\@tsp{3} % triple space
}
\DeclareOptionX{draftspace}{
  \setlength\overfullrule{5pt}
  \def\DUT@toc@sep{0} % content page separation (# blanklines between entries)
  \def\@ssp{1} % single space
  \def\@dsp{1} % double space
  \def\@tsp{2} % triple space
}
\DeclareOptionX{pagebreak}{
  \@DUT@pagebreaktrue
  \PassOptionsToClass{openright}{book}
}
\DeclareOptionX{nopagebreak}{
  \@DUT@pagebreakfalse
  \PassOptionsToClass{openany}{book}
}
\DeclareOptionX{dottedpart}{\@dottedparttrue}
\DeclareOptionX{dottedchapter}{\@dottedchaptertrue}
\DeclareOptionX{boldpart}{\@boldparttrue}
\DeclareOptionX{boldchapter}{\@boldchaptertrue}
\DeclareOptionX{dottedtoc}{
  \setkeys{drexel-thesis.cls}{dottedpart}
  \setkeys{drexel-thesis.cls}{dottedchapter}
}
\DeclareOptionX{boldtoc}{
  \setkeys{drexel-thesis.cls}{boldpart}
  \setkeys{drexel-thesis.cls}{boldchapter}
}
\DeclareOptionX{indenttoc}{\@indenttoctrue}
\DeclareOptionX{blackref}{\@blackreftrue}
\DeclareOptionX{draftwatermark}{%
  \@draftwatermarktrue
  \@xdraftwatermarkfalse}
\DeclareOptionX{xdraftwatermark}{%
  \@draftwatermarkfalse
  \@xdraftwatermarktrue}
 \DeclareOptionX{nodraftwatermark}{%
  \@draftwatermarkfalse
  \@xdraftwatermarkfalse}
\DeclareOptionX{indentfirst}{\@indentfirsttrue}
\DeclareOptionX{oneside}{\PassOptionsToClass{oneside}{book}}
\DeclareOptionX{twoside}{\PassOptionsToClass{twoside}{book}}
\DeclareOptionX{final}{
  \@finaltrue
  \setkeys{drexel-thesis.cls}{oneside}
  \setkeys{drexel-thesis.cls}{finalspace}
  \setkeys{drexel-thesis.cls}{pagebreak}
  \setkeys{drexel-thesis.cls}{dottedtoc}
}
\DeclareOptionX{draft}{
  \@finalfalse
  \setkeys{drexel-thesis.cls}{twoside}
  \setkeys{drexel-thesis.cls}{draftspace}
  \setkeys{drexel-thesis.cls}{nopagebreak}
  \setkeys{drexel-thesis.cls}{dottedtoc}
  \setkeys{drexel-thesis.cls}{xdraftwatermark}
}
\DeclareOptionX{daring}{
  \@daringtrue
}
\DeclareOptionX{approvalform}{
  \@approvalformtrue
}
\DeclareOptionX*{\PassOptionsToClass{\CurrentOption}{book}}
\ExecuteOptionsX{final}
\ProcessOptionsX\relax
\LoadClass{book}
\RequirePackage{setspace}
\RequirePackage{calc}
\RequirePackage[final]{graphicx}
\if@draftwatermark
  \RequirePackage[draft,allpages,mark={DRAFT \today},
                  fontsize=2cm,grayness=0.9]{draftmark}
\fi
\if@xdraftwatermark
  \RequirePackage{eso-pic}
  \RequirePackage{xcolor}
  \definecolor{watermark}{rgb}{0.9,0.9,0.9}
  \newcommand{\watermarkfont}{%
    \color{watermark}}%\fontfamily{phv}\selectfont\sffamily}
  \newcommand{\watermark}{%
    \rule{4.25in}{0pt}\raisebox{5.5in}{%
      \makebox[0pt]{%
\rotatebox[origin=c]{45}{%
          \resizebox{!}{1.5in}{%
            {\watermarkfont DRAFT}%
    }}}}
    \raisebox{4.2in}{%
      \makebox[0pt]{%
\rotatebox[origin=c]{45}{%
          \resizebox{!}{.6in}{%
            {\watermarkfont \qquad\today}%
    }}}}}
  \AddToShipoutPicture{\watermark{}}
\fi
\if@indentfirst
  \RequirePackage{indentfirst}
\fi
\RequirePackage[letterpaper,margin=1in,bindingoffset=0.5in]{geometry}
\if@DUT@subfigopt
  \RequirePackage[subfigure]{tocloft}
\else
  \RequirePackage{tocloft}
\fi
\RequirePackage[margin=10pt,labelfont=bf,labelsep=colon]{caption}
\RequirePackage{floatrow}
\floatsetup[table]{style=plaintop}
\floatsetup{captionskip=10pt}
\floatsetup[table]{captionskip=10pt}
\if@DUT@subfigopt
  \RequirePackage{fr-subfig}
\fi
\if@DUT@subfigopt
  \RequirePackage{subfig}
\fi
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}
% \renewcommand{\footrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0pt}
% \fancyfoot[RE,LO]{\scshape\leftmark}
\fancyfoot[CE,CO]{\thepage}
% \fancyfoot[LE,RO]{\scshape\rightmark}
% \fancyhead[LE,RO]{\thepage}
\fancypagestyle{plain}{
  \fancyhf{}
  \renewcommand{\footrulewidth}{0pt}
  % \fancyhead[LE,RO]{\thepage}
  \fancyfoot[CE,CO]{\thepage}

}
\fancypagestyle{empty}{
  \fancyhf{}
  \renewcommand{\footrulewidth}{0pt}
}
\RequirePackage{ifpdf}
\if@approvalform\RequirePackage{tikz}\fi
\if@approvalform\RequirePackage{pdfpages}\fi
\renewcommand{\title}[1]{
  \gdef\@DUT@title{#1}
  \ifpdf \hypersetup{pdftitle={\@DUT@title}} \fi
}
\let\DUTtitle\title
\def\@DUT@title{\@latex@error{No \noexpand\title given}\@ehc}
\renewcommand{\author}[1]{
  \gdef\@DUT@author{#1}
  \ifpdf \hypersetup{pdfauthor={\protect\@DUT@author}} \fi
}
\let\DUTauthor\author
\def\@DUT@author{\@latex@error{No \noexpand\author given}\@ehc}

%adding option to place approval form location
\newcommand{\approvalForm}[1]{\def\@approvalForm{#1}}
\approvalForm{docs/approval-form.pdf}
\def\DUTmonth#1{\gdef\@DUT@month{#1}}
\def\@DUT@month{\@latex@error{No \noexpand\DUTmonth given}\@ehc}
\def\DUTyear#1{\gdef\@DUT@year{#1}}
\def\@DUT@year{\@latex@error{No \noexpand\DUTyear given}\@ehc}
\def\degree#1{\gdef\@DUT@degree{#1}
  \ifpdf \hypersetup{pdfsubject={\@DUT@degree}} \fi}
\let\DUTdegree\degree
\def\@DUT@degree{\@latex@error{No \noexpand\degree given}\@ehc}
\def\copyrighttext#1{\gdef\@DUT@copyrighttext{#1}}
\let\DUTcopyrighttext\copyrighttext
\def\@DUT@copyrighttext{All Rights Reserved.}
\newcounter{@DUT@advisors}
\setcounter{@DUT@advisors}{0}
\newcommand\advisor[1]{%
  \ifx
    \empty{#1}
  \else
    \stepcounter{@DUT@advisors}
    \expandafter\gdef\csname @DUT@advisor@\alph{@DUT@advisors}\endcsname{#1}
  \fi}
\newcommand\@DUT@advisor{%
  \ifnum\c@@DUT@advisors=0%
    \@latex@error{No \noexpand\advisor given}\@ehc%
  \else%
    \begingroup%
      \count@\@ne
      \csname @DUT@advisor@\@alph\count@\endcsname\relax%
      \loop\ifnum\count@<\c@@DUT@advisors%
        \advance\count@\@ne
        \ifnum\count@=\c@@DUT@advisors%
          \ifnum\c@@DUT@advisors=2%
            {\ and\ }%
          \else%
            {, and\ }%
          \fi%
        \else%
          {,\ }%
        \fi%
        \csname @DUT@advisor@\@alph\count@\endcsname\relax%
      \repeat%
    \endgroup%
  \fi}
\newlength\@DUT@baselineskip
\setlength\@DUT@baselineskip\baselineskip
\newcommand{\set@lines@length}[2]{%
  \setlength{#1}{\@DUT@baselineskip * {#2}}}
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
    \let\@svsec\@empty
  \else
    \refstepcounter{#1}%
    \protected@edef\@svsec{\@seccntformat{#1}\relax}%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      \setstretch{\@ssp}%
      #6{%
        \@hangfrom{\hskip #3\relax\@svsec}%
          \interlinepenalty \@M #8\@@par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth \else
        \protect\numberline{\csname the#1\endcsname}%
      \fi
      #7}%
  \else
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec #8}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}
\def\@ssect#1#2#3#4#5{%
  \@tempskipa #3\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      \setstretch{\@ssp}%
      #4{%
        \@hangfrom{\hskip #1}%
          \interlinepenalty \@M #5\@@par}%
    \endgroup
  \else
    \def\@svsechd{#4{\hskip #1\relax #5}}%
  \fi
  \@xsect{#3}}
\newcommand\@DUT@clearpage{%
  \if@DUT@pagebreak \clearpage \else \blanklinesp{4} \fi}
\newcommand\@DUT@cleardoublepage{%
  \if@DUT@pagebreak \cleardoublepage \else \blanklinesp{4} \fi}
\newcommand\@DUT@newpage{%
  \if@DUT@pagebreak \newpage \else \blanklinesp{4} \fi}
\newcommand\@DUT@blankpage{%
  \if@DUT@pagebreak
    \@DUT@clearpage
    \thispagestyle{empty}
    \hphantom{0}
    \@DUT@clearpage
  \fi
}
\newlength\@DUT@blanklines@length
\newcommand{\@blanklines@setup}[1]{%
  \set@lines@length{\@DUT@blanklines@length}{#1}%
  %\typeout{BLANKLINES -- #1}%
  %\showthe\@DUT@blanklines@length
}
\newcommand{\blanklines}{\@ifstar\@sblanklines\@blanklines}
\newcommand{\@blanklines}[1]{%
  \@blanklines@setup{#1}%
  \vspace{\@DUT@blanklines@length}}
\newcommand{\@sblanklines}[1]{%
  \@blanklines@setup{#1}%
  \vspace*{\@DUT@blanklines@length}}
\newcommand{\blanklinesp}{\@ifstar\@sblanklines\@blanklines}
\newcommand{\@blanklinesp}[1]{%
  \@blanklines@setup{#1}%
  \vspace{\@DUT@blanklines@length}\par}
\newcommand{\@sblanklinesp}[1]{%
  \@blanklines@setup{#1}%
  \vspace*{\@DUT@blanklines@length}\par}
\RequirePackage[final]{hyperref}
\hypersetup{colorlinks}
\if@blackref
  \hypersetup{linkcolor=black}
  \hypersetup{anchorcolor=black}
  \hypersetup{citecolor=black}
  \hypersetup{filecolor=black}
  \hypersetup{menucolor=black}
  \hypersetup{runcolor=black}
  \hypersetup{urlcolor=black}
\fi
\RequirePackage[all]{hypcap}
\set@lines@length{\floatsep}{\@ssp}
\set@lines@length{\textfloatsep}{\@tsp}
\newcommand\@DUT@maketitle{
  \begin{@DUT@titlepage}%
    \ifpdf
      \pdfbookmark[-1]{Front Matter}{Front Matter}
      \pdfbookmark[0]{Title Page}{Title Page}
    \fi
    \if@DUT@pagebreak \vspace*{\fill} \fi
    \begin{center}%
      {\setstretch{\@tsp} \bfseries \@DUT@title \par}%
      \blanklines{\@tsp}%
      \setstretch{\@ssp}%
      A Dissertation submitted in partial fulfillment of the requirements for the degree of 
      Doctor of Philosophy at George Mason University\\
      \blanklines{\@tsp}%
      \setstretch{\@ssp}%
      by\\
      \blanklines{\@tsp}%
      \setstretch{\@ssp}%
      \@DUT@author\\
      Master of Science\\
      Georgetown University, 2012\\
      Bachelor of Arts\\
      The College of William and Mary, 2009\\
      \blanklines{\@tsp}%
      \setstretch{\@ssp}%
      Director: Jason Kinser, Chair\\ 
      Department of Computational and Data Sciences\\
      \blanklines{\@tsp}%
      \setstretch{\@ssp}%
      \@DUT@month\ \@DUT@year\\
      George Mason University\\
      Fairfax, Virginia\\
    \end{center}\par
    \if@DUT@pagebreak \vspace*{\fill}
    \centering
    \mbox{\smash{\includegraphics[width=4cm]{figure/drexel-logo}}}
    \fi
  \end{@DUT@titlepage}%
  \setcounter{footnote}{0}%
  \global\let\maketitle\relax
  \global\let\and\relax
}
\newcommand\@DUT@makecopyright{%
  \@DUT@clearpage
  \ifpdf \pdfbookmark[0]{Copyright Page}{Copyright Page} \fi
  \if@DUT@pagebreak \vspace*{\fill} \fi
  \begin{center}
    \copyright~Copyright \@DUT@year \hspace{1mm} \@DUT@author.\\
    \@DUT@copyrighttext
    \setstretch{\@ssp}%
  \end{center}
  \if@DUT@pagebreak \vspace*{\fill} \fi
  \thispagestyle{empty}
}
\newcommand\copyrighttextCCBYSA{
  \begin{quote}
    This work is licensed under the terms of the Creative Commons
    Attribution-ShareAlike 4.0 International license.  The license is
    available at \\
    \url{http://creativecommons.org/licenses/by-sa/4.0/}.
  \end{quote}
}
\let\DUTcopyrighttextCCBYSA\copyrighttextCCBYSA
\renewcommand\chaptermark[1]{%
  \markboth {%
    \ifnum \c@secnumdepth >\m@ne
      \if@mainmatter
        \@chapapp\ \thechapter: \ %
      \fi
    \fi
    #1}{}}
\newcommand\schaptermark[1]{\markboth{#1}{}}
\renewcommand{\sectionmark}[1]{%
  \markright {%
    \ifnum \c@secnumdepth >\z@
      \thesection \ %
    \fi
    #1}}
\renewcommand\mainmatter{%
  \@DUT@cleardoublepage
  \@mainmattertrue
  \pagenumbering{arabic}}
\renewcommand\part{%
  \if@openright
    \@DUT@cleardoublepage
  \else
    \@DUT@clearpage
  \fi
  \thispagestyle{plain}%
  \@tempswafalse
  \if@DUT@pagebreak \null\vfil \fi
  \secdef\@part\@spart}
\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >-2\relax
      \refstepcounter{part}%
      \addcontentsline{toc}{part}%
        {\protect\numberline{\thepart}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    \markboth{}{}%
    {\centering
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >-2\relax
       \large\bfseries \partname\nobreakspace\thepart
       \par
       \vskip 20\p@
     \fi
     \large \bfseries #2\par}%
    \@endpart}
\def\@spart#1{%
    {\centering
     \interlinepenalty \@M
     \normalfont
     \large \bfseries #1\par}%
    \@endpart}
\def\@endpart{
  \if@DUT@pagebreak \vfil\newpage \fi
  \if@twoside
    \if@openright
      \null
      \thispagestyle{empty}%
      \newpage
    \fi
  \fi}
\renewcommand\chapter{%
  \if@openright\@DUT@cleardoublepage\else\@DUT@clearpage\fi
  \thispagestyle{plain}%
  \global\@topnum\z@
  \@afterindentfalse
  \secdef\@chapter\@schapter}
\newcommand\listed@schapter[1]{%
  \if@openright\@DUT@cleardoublepage\else\@DUT@clearpage\fi
  \thispagestyle{plain}%
  \global\@topnum\z@
  \phantomsection
  \addcontentsline{toc}{chapter}{#1}
  \@afterindentfalse
  \@schapter{#1}}
\def\@chapter[#1]#2{
  \let\@DUT@label=\label%
  \renewcommand{\label}[1]{}%
  \typeout{\@chapapp\space\thechapter. #1}%
  \let\label=\@DUT@label%
  \ifnum \c@secnumdepth >\m@ne
    \if@mainmatter
      \refstepcounter{chapter}%
      \addcontentsline{toc}{chapter}{\protect\numberline{\thechapter}#1}%
    \else
      \addcontentsline{toc}{chapter}{#1}%
    \fi
  \else
    \addcontentsline{toc}{chapter}{#1}%
  \fi
  \chaptermark{#1}%
  \@makechapterhead{#2}%
  \@afterheading}
\def\@makechapterhead#1{%
  \blanklines*{-1}% HACK!
  {\parindent \z@ \raggedright \normalfont
    \setstretch{\@ssp}%
    \interlinepenalty\@M
    \ifnum \c@secnumdepth >\m@ne
      \if@mainmatter
        \center \large \bfseries \@chapapp{} \thechapter: #1\par\nobreak
      \else
        \center \large \bfseries #1\par\nobreak
      \fi
    \else
      \center \large \bfseries #1\par\nobreak
    \fi
    \blanklines{4}
  }
}
\def\@schapter#1{
  \typeout{\@chapapp\space #1}%
  \schaptermark{#1}%
  \@makeschapterhead{#1}%
  \@afterheading}
\def\@makeschapterhead#1{%
  \blanklines*{-1}% HACK!
  {\parindent \z@ \raggedright \normalfont
    \setstretch{\@ssp}%
    \interlinepenalty\@M
    \center \large \bfseries  #1\par\nobreak
    \blanklines{4}%
  }
}
\renewcommand\section{%
  \@startsection {section}{1}{\z@}%
    {-3.5ex \@plus -1ex \@minus -.2ex}%
    {2.3ex \@plus.2ex}%
    {\normalfont\large\bfseries}}
\newenvironment{preamble}{%
  \cleardoublepage
  \@mainmatterfalse
  \pagenumbering{roman}
  \if@approvalform
    \ifpdf
      \pdfbookmark[-1]{Approval Form}{Approval Form}
    \fi
    \DUT@makeapprovalform
  \fi
  \@DUT@maketitle
  \@DUT@makecopyright
}{%
  \@DUT@blankpage
}
\let\DUTpreamble\preamble
\let\endDUTpreamble\endpreamble
\newenvironment{dedications}{%
  \@DUT@clearpage
  \ifpdf
    \phantomsection
    \pdfbookmark[0]{\DUTdedicationsname}{\DUTdedicationsname}
  \fi
  \chapter*{\DUTdedicationsname}
  \setstretch{\@ssp}%
}{%
  \par\setstretch{\@ssp}%
}
\let\DUTdedications\dedications
\let\endDUTdedications\enddedications
\newenvironment{acknowledgments}{%
  \@DUT@clearpage
  \ifpdf
    \phantomsection
    \pdfbookmark[0]{\DUTacknowledgmentsname}{\DUTacknowledgmentsname}
  \fi
  \chapter*{\acknowledgmentsname}
% \setstretch{\@dsp}%
\setstretch{\@ssp}%
}{%
  \par\setstretch{\@ssp}%
}
\let\DUTacknowledgments\acknowledgments
\let\endDUTacknowledgments\endacknowledgments
\newenvironment{abstract}{%
  \listed@schapter{\abstractname}%
  \blanklines{-4}%
    % \begin{center}
      \setstretch{\@ssp}%
      \blanklines{3}%
      \MakeUppercase{\@DUT@title}\\
      \@DUT@author
      , PhD.\\
      George Mason University, 2023\\
      \ifdaring{%
        \ifnum\c@@DUT@advisors=\@ne%
        Advisor:
        \else%
        Advisors:
        \fi}{}
      \@DUT@advisor\\
    % \end{center}
  % \blanklines{\@tsp}%
  % \setstretch{\@dsp}%
  \blanklines{4}%
  \setstretch{\@dsp}%
  \@nobreaktrue
  \@afterindentfalse
  \@afterheading
}{%
  \par\setstretch{\@ssp}%
}
\let\DUTabstract\abstract
\let\endDUTabstract\endabstract
\newenvironment{thesis}{%
  \mainmatter
  \setstretch{\@dsp}%
}{%
  \setstretch{\@ssp}%
  \@DUT@clearpage
  \ifpdf
    \phantomsection
    \pdfbookmark[-1]{Back Matter}{Back Matter}
  \fi
}
\let\DUTthesis\thesis
\let\endDUTthesis\endthesis
\newenvironment{vita}{%
  \listed@schapter{\vitaname} % commented out to take out of TOC
  % \setstretch{\@dsp}% % does this make pages double spaced?
  \setstretch{\@ssp}
}{%
  \par\setstretch{\@ssp}%
}
\let\DUTvita\vita
\let\endDUTvita\endvita
\renewenvironment{titlepage}%
  {\thispagestyle{empty}\setcounter{page}\z@}%
  {\@DUT@newpage}
\let\@DUT@titlepage\titlepage
\let\end@DUT@titlepage\endtitlepage
\let\@DUT@appendix\appendix
\renewcommand\appendix{%
  \@DUT@appendix%
  \addtocontents{toc}{%
    \protect\cftsetindents{chapter}{0pt}{6.5em}
    \protect\renewcommand{\protect\cftchappresnum}{\protect\appendixname~}
    \protect\renewcommand{\protect\cftchapaftersnum}{:}}
  \setstretch{\@dsp}%
}
\let\DUTappendix\appendix
\renewcommand\cfttoctitlefont{%
  \@DUT@clearpage
  \ifpdf
    \phantomsection
    \pdfbookmark[0]{\contentsname}{\contentsname}
  \fi
  \chapter*}
\setlength\cftaftertoctitleskip{0pt}
\newlength\l@DUT@toc@sep
\set@lines@length{\l@DUT@toc@sep}{\@DUT@toc@sep}
\newcommand{\DUT@toc@high@font}{\normalsize\scshape}
\newcommand{\DUT@toc@low@font}{\normalfont}
\newcommand{\DUT@toc@page@font}{\normalfont}
\setlength\cftbeforepartskip{\l@DUT@toc@sep}
\cftsetindents{part}{0pt}{3em}
\if@boldpart   \renewcommand{\cftpartfont}{\bfseries\DUT@toc@high@font}
\else          \renewcommand{\cftpartfont}{\DUT@toc@high@font}
\renewcommand{\cftpartleader}{\cftdotfill{\cftpartdotsep}}
\if@dottedpart \renewcommand{\cftpartdotsep}{\cftdotsep}
\else          \renewcommand{\cftpartdotsep}{\cftnodots} \fi
\renewcommand{\cftpartpagefont}{\DUT@toc@page@font}
\setlength\cftbeforechapskip{\l@DUT@toc@sep}
\cftsetindents{chapter}{0pt}{1.5em}
\renewcommand{\cftchapaftersnum}{.}
\if@boldchapter   \renewcommand{\cftchapfont}{\bfseries\DUT@toc@high@font}
\else             \renewcommand{\cftchapfont}{\DUT@toc@high@font} \fi
\renewcommand{\cftchapleader}{\cftdotfill{\cftchapdotsep}}
\if@dottedchapter \renewcommand{\cftchapdotsep}{\cftdotsep}
\else             \renewcommand{\cftchapdotsep}{\cftnodots} \fi
\renewcommand{\cftchappagefont}{\DUT@toc@page@font}
\setlength\cftbeforesecskip{\l@DUT@toc@sep}
\cftsetindents{section}{0pt}{3em}
\renewcommand{\cftsecfont}{\DUT@toc@low@font}
\renewcommand{\cftsecpagefont}{\DUT@toc@page@font}
\setlength\cftbeforesubsecskip{\l@DUT@toc@sep}
\cftsetindents{subsection}{0pt}{4.5em}
\renewcommand{\cftsubsecfont}{\DUT@toc@low@font}
\renewcommand{\cftsubsecpagefont}{\DUT@toc@page@font}
\setlength\cftbeforesubsubsecskip{\l@DUT@toc@sep}
\cftsetindents{subsubsection}{0pt}{6em}
\renewcommand{\cftsubsubsecfont}{\DUT@toc@low@font}
\renewcommand{\cftsubsubsecpagefont}{\DUT@toc@page@font}
\if@indenttoc
  \cftsetindents{chapter}{0pt}{1.5em}
  \cftsetindents{section}{1.5em}{2.3em}
  \cftsetindents{subsection}{3.8em}{3.2em}
  \cftsetindents{subsubsection}{7.0em}{4.1em}
\fi
\renewcommand\cftloftitlefont{\listed@schapter}
\setlength\cftafterloftitleskip{-12pt}
\setlength{\cftbeforefigskip}{\l@DUT@toc@sep}
\setlength{\cftfigindent}{0pt}
\renewcommand{\cftfigfont}{\DUT@toc@low@font}
\renewcommand{\cftfigpagefont}{\DUT@toc@page@font}
\renewcommand\cftlottitlefont{\listed@schapter}
\setlength\cftafterlottitleskip{-12pt}
\setlength{\cfttabindent}{0pt}
\setlength{\cftbeforetabskip}{\l@DUT@toc@sep}
\renewcommand{\cfttabfont}{\DUT@toc@low@font}
\renewcommand{\cfttabpagefont}{\DUT@toc@page@font}
\renewenvironment{thebibliography}[1]
     {\bib@heading%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand{\bib@heading}{%
  \listed@schapter{\bibname}%
}
\renewenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
  \listed@schapter{\indexname}
  \global\columnwidth\textwidth
  \global\advance\columnwidth-\columnsep
  \global\divide\columnwidth\tw@
  \global\hsize\columnwidth
  \global\linewidth\columnwidth
  \global\@twocolumntrue
  \global\@firstcolumntrue
  \col@number \tw@
                \parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \columnseprule \z@
                \columnsep 35\p@
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\@DUT@clearpage\fi}
\newcommand\dedicationsname{Dedications}
\let\DUTdedicationsname\dedicationsname
\newcommand\acknowledgmentsname{Acknowledgments}
\let\DUTacknowledgmentsname\acknowledgmentsname
\renewcommand\contentsname{Table of Contents}
\let\DUTcontentsname\contentsname
\renewcommand\listtablename{List of Tables}
\let\DUTlisttablename\listtablename
\renewcommand\listfigurename{List of Figures}
\let\DUTlistfigurename\listfigurename
\newcommand\abstractname{Abstract}
\let\DUTabstractname\abstractname
\renewcommand\bibname{Bibliography}
\let\DUTbibname\bibname
\renewcommand\appendixname{Appendix}
\let\DUTappendixname\appendixname
\newcommand\vitaname{Biography}
\let\DUTvitaname\vitaname
\newcommand{\DUTnames}{%
  \gdef\dedicationsname{\DUTdedicationsname}%
  \gdef\contentsname{\DUTcontentsname}%
  \gdef\listtablename{\DUTlisttablename}%
  \gdef\listfigurename{\DUTlistfigurename}%
  \gdef\bibname{\DUTbibname}%
  \gdef\appendixname{\DUTappendixname}%
  \gdef\vitaname{\DUTvitaname}%
}
\AtEndDocument{\@DUT@blankpage}
\endinput
%%
%% End of file `drexel-thesis.cls'.
