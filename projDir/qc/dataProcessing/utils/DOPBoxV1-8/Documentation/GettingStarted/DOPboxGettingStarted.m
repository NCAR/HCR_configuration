%% Header to define the Title and authors.
% \documentclass[12pt]{article}
%
% \title{\textbf{DOPbox}\\ Discrete Orthogonal Polynomial Toolbox\\
% Getting Started}
% 
% \author{Matthew Harker and Paul O'Leary\\
% Institute for Automation\\
% University of Leoben\\
% A-8700 Leoben,
% Austria\\
% URL: automation.unileoben.ac.at\\
% \\
% Original: January 9, 2013\\
% $\copyright$ 2013\\
% \\
% Last Modified: \today}
%
%%
%\section{Toolbox description} 
%
% This package contains a set of m-files which implement a
% \textbf{D}iscrete \textbf{O}rthogonal \textbf{P}olynomial toolbox
% called
% \textbf{DOPbox}. There are many applications of discrete
% opthogonal polynomials, e.g., in the solution of ordinary
% differential equations (ODEs), partial differential equations
% (PDEs), boundary value problems (BVPs), initial value problems
% (IVPs), inverse-BVPs etc.
%
% This documentation describes the functions briefly and gives examples as
% to how they are used. The theory behind this toolbox can be found in the
% following publications: The paper~\cite{oleary2008b}
%{\footnotesize
% \begin{verbatim}
% @inproceedings{
% oleary2008b,
%    Author = {O'Leary, P. and Harker, M.},
%    Title = {An Algebraic Framework for Discrete Basis Functions in Computer Vision},
%    BookTitle = {2008 $6^{\textrm{th}}$ ICVGIP},
%    Address= {Bhubaneswar, India},
%    Publisher = {IEEE},
%    Pages = {150-157},
%  Year = {2008} }
% \end{verbatim}
% }
%
% provides details of the synthesis algorithm and the application of
% the Gram polynomials to image processing. The paper~\cite{oleary2010C}
%
% {\footnotesize
% \begin{verbatim}
% @inproceedings{oleary2010C,
%    Author = {O'Leary, P. and Harker, M.},
%    Title = {Discrete Polynomial Moments and Savitzky-Golay Smoothing},
%    BookTitle = {Waset Special Journal},
%    Volume = {72},
%    DOI = {},
%    Pages = {439--443},
% Year = {2010}}
% \end{verbatim}
% }
%
% provides information on when the Gram-Schmidt orthogonalization fails and
% also the implementation of local polynomial approximations. The
% application of \lstinline{DOPbox} to BVPs and IVBPs is documented in~\cite{Oleary2012},
%
%{\footnotesize
% \begin{verbatim}
% @article{Oleary2012,
%   author    = {Paul O'Leary and
%                Matthew Harker},
%   title     = {A Framework for the Evaluation of Inclinometer Data in the
%                Measurement of Structures},
%   journal   = {IEEE T. Instrumentation and Measurement},
%   volume    = {61},
%   number    = {5},
%   year      = {2012},
%   pages     = {1237-1251},
%   ee        = {http://dx.doi.org/10.1109/TIM.2011.2180969}}
% \end{verbatim}
% }
%
%%
%\section{Toolbox Functions}
%
% All the m-files contain headers with help for the specific function,
% see for example:
help dopNodes
%%
% Clean up prior to running the code 
%
clear all;
close all;
%%======================================================================
%\subsection{Generate a set of nodes -- \lstinline{dopNodes.m}}
% this function generates a set of nodes which can be used to synthesize a
% set of basis functions, see \lstinline{help dopNodes} for details of the
% types of nodes supported.
%
% The function is used in the following manner
%
noNodes = 10;
%
% generate the Gram and Chebyshev nodes
%
gNodes = dopNodes( noNodes, 'Gram' );
cNodes = dopNodes( noNodes, 'Cheby' );
%
% Generate the modiied nodes
%
geNodes = dopNodes( noNodes, 'GramEnds' );
ceNodes = dopNodes( noNodes, 'ChebyEnds' );
%%
% It is important to note that the nodes need to be generated exactly at
% the points where the constraints are to be applied, in particular for
% boundary values problems. Both the Gram and Chebyshev nodes do not reach
% to the end of the support; this would  mean that constraints could
% not be placed exactly at the ends of the support. For this reason
% modified Gram and Chenyshev nodes are also provided, they are
% modified so that they start at $-1$ and end at $1$.
%
%
fig1 = figure;
a(1) = subplot(4,1,1);
plot( gNodes, zeros(size(gNodes)), 'ko', 'MarkerFaceColor', 'k');
hold on;
plot([-1,1],[0,0], 'k');
axis([-1,1,-1,1]);
ylabel('Gram ');
%
a(3) = subplot(4,1,2);
plot( cNodes, zeros(size(gNodes)), 'ko', 'MarkerFaceColor', 'k');
hold on;
plot([-1,1],[0,0], 'k');
axis([-1,1,-1,1]);
ylabel('Cheby');
%
a(2) = subplot(4,1,3);
plot( geNodes, zeros(size(gNodes)), 'ko', 'MarkerFaceColor', 'k');
hold on;
plot([-1,1],[0,0], 'k');
axis([-1,1,-1,1]);
ylabel('GramE');
%
a(4) = subplot(4,1,4);
plot( ceNodes, zeros(size(gNodes)), 'ko', 'MarkerFaceColor', 'k');
hold on;
plot([-1,1],[0,0], 'k');
axis([-1,1,-1,1]);
%
xlabel( 'Support' );
ylabel('ChebyE');
%
linkaxes( a, 'x');
%
%\caption{Node placements generated by \lstinline{dopNodes} with the options
% for the Gram, Cheby, GramEnds and ChebyEnds.}
%%
% In the next figure the end of the support is shown.
%
figure(fig1);
axis([-1,-0.65,-1,1]);
%
%\caption{Node placements at the end of the support as
% generated by \lstinline{dopNodes} with the options
% for the Gram, Cheby, GramEnds and ChebyEnds.}
%%======================================================================
% \subsection{Generate a set of discrete orthogonal polynomials -- \lstinline{dop.m}}
%
% The function \lstinline{dop.m} generates a set of discrete orthogonal
% polynomials. It has multiple options (see \lstinline{help dop} for
% details), the most important features are: the fuinction generates a set
% of ortho-normal basis functions from a set of nodes. It can also
% optionally generate the differentials of the basis functions and the
% coefficients for later interpolation with the basis functions.
%
% The generation of a set of Gram polynomials from the Gram nodes is
% now demonstrated. A large number of nodes are used to show the nature of
% the basis functions.
%
noNodes = 100;
noBfs = noNodes;
%
% generate a complete set of basis functions.
%
gNodes = dopNodes( noNodes, 'Gram' );
B = dop( gNodes, noBfs );
%
fig2 = figure;
for k=1:5
    plot( gNodes, B(:,k), 'k' );
    hold on;
end;
axis([-1,1,-0.3, 0.3]);
grid on;
xlabel('Support');
ylabel('Value');
%
% \caption{The first $5$ Gram basis functions.}
%%
% \subsubsection{Properties of the Basis functions}
%
% The basis functions $\M{B}$ generated by this function are ortho-normal, i.e.
% the Gram matrix for the basis functions
% \begin{equation}
%   \M{G} \defas \MT{B} \, \M{B} = \M{I}
% \end{equation}
% should be the identity matrix. Consequently, the projection onto the 
% orthogonal complement of $\M{G}^\perp = \M{G} - \M{I} = \M{0}$ should be
% exactly a matrix of zeros.
% However, due to numerical errors this is not the case.  $\M{G}^\perp$ can
% be inspected to see the quality of the basis functions.
%
G = B' * B;
Gperp = G - eye( noBfs );
%
% Now visualizing the matrix Gperp,
%
xScale = (1:noBfs) - 1;
yScale = xScale;
fig3 = figure;
imagesc( xScale, yScale, Gperp );
colorbar;
xlabel('Degree $$d$$');
ylabel('Degree $$d$$');
axis image;
% \caption{The projection onto the orthogonal complement of the Gram matrix, 
% computed for a set of $100$ Gram basis functions. The matrix should be exactly zeror.}
%%
%
% Ideally the matrix $\M{G}^\perp$ should be exactly zero. The Frobenius
% norm $\| \M{G} \|_f$, i.e. the sum of the squares of all entries, can be used as an
% estimate for the quality of the basis functions.
%
fNormGperp = norm( Gperp, 'fro')
%%
% The number of significant digits can be estimated as $n_s = - \log_{10}(\| \M{G} \|_f)$
%
ns = - log10( fNormGperp )

%%
% As can be seen the error is very small. Consequently, the quality of the
% basis functions is very high.
%
%%
% \subsubsection{Generating complex basis functions}
%
% Consider the case of wishing to generate a set of basis functions from an
% arbitrary set of nodes in the complex plane. This, for example, may be
% desirable when wishing to modell an damped oscillator.
%
% In the example presented here we wish to generate a set of complex basis
% functions with a predefined envelope. The complex basis functions should
% be in quadrature to each other, this insures local shift invariance.
%
load arbitraryData;
%
n = length(r);
x = linspace(0,1,n)';
%
figure;
plot( x, r, 'k');
grid on;
xlabel('Support');
ylabel('Envelope');
%\caption{Envelope for the desired basis functions}
%%
% Now defining the number of cycles required for the basis functions
%
noCycles = 10;
%
% and compute the corresponding angles
%
phi = linspace(0,2*pi*noCycles,n )';
%
% Compute the rean and imaginary components of the node placements
%
nr = r.*cos( phi );
ni = r.*sin( phi );
n = nr + i * ni;
%
figure;
plot( nr, ni, 'k');
hold on;
plot( nr, ni, 'k.');
grid on;
xlabel('Real');
ylabel('Imaginary');
axis equal;
%\caption{Node placements for the synthesis of the complex basis functions.}
%%
% Now synthesize the basis functions (only the first 3 are computed for
% simplicity here.
%
noBfs = 3;
B = dop( n, noBfs );
%
figure;
subplot(2,1,1);
plot( x, real(B(:,2)) , 'r');
hold on;
plot( x, imag(B(:,2)) , 'b');
grid on
ylabel( 'B(:,2)');
plot(x,r/norm(r),'k');
plot(x,-r/norm(r),'k');
%
subplot(2,1,2);
plot( x, real(B(:,3)) , 'r');
hold on;
plot( x, imag(B(:,3)) , 'b');
grid on
ylabel( 'B(:,3)');
%
xlabel('Support');
%\caption{Sets of complex basis functions synthesized using \lstinline{dop}}

%%======================================================================
% \subsubsection{Using \lstinline{dop.m} to compute a regularizing differential operator.}
% 
% In this example a polynomial is used, since the function and its
% derivatives can be computed analytically. This permits an objective
% comparison of the numerically computed and the analytical derivatives.
%

% Generate a polynomial function for test purposes:
noPts = 100;
x = linspace( -1, 1, noPts )';
p = [5,0,-5,0,-5,0];
dp = polyder( p );
y = polyval( p, x );
dy = polyval( dp, x );
%
% Plot the analytical polynomials and their derivatives together with the
% numerically computed values
%
fig4 = figure;
subplot(2,1,1);
plot( x, y, 'k');
ylabel( '$$y(x)$$' );
grid on;
%
subplot(2,1,2);
plot( x, dy, 'k');
ylabel( '$$d y(x)/dx$$' );
grid on;
%
% \caption{the test function and its analytical derivative}
%%
%
% Now investigating the properties of numerical approximations to
% derivatives: the MATLAB matrix implicitly used to numerically compute
% derivatives can be extracted as follows:
%
[Dx, Dy] = gradient( eye(noPts));
%
% correct for the step size
%
h = x(2) - x(1);
DyMat = Dy / h ;
%
% and using this matrix to estimate the derivative of y(x).
%
dyMAT = DyMat * y ;
%%
% The function dop.m can also generate a derivative of the basis functions.
% Given the basis functions $\M{B}$ and their derivatives $\dot{\M{B}}$ a
% differential operator $\M{D}$ can be computed:
% \begin{equation}
%   \dot{\M{B}} = \M{D} \, \M{B}
% \end{equation}
% consequently
% \begin{equation}
%   \M{D}_r \defas \M{D} \, \M{B} \,\MT{B}  = \dot{\M{B}} \, \MT{B}.
% \end{equation}
%
% This is a regularizing differential operator, it is the differential
% operator $\M{D}$ projected onto the basis functions $\M{B}$
%
[B, dB] = dop( x, length(p) );
Dr = dB * B';
%
% and using D to estimate the derivatives
%
dyDop = Dr * y;
%
fig5 = figure;
plot( x, dy - dyMAT, 'k');
hold on;
plot( x, dy - dyDop, 'r');
grid on;
xlabel('Support');
ylabel('Error in $$d y(x)$$');
legend('Matlab Gradient','dop derivative','location','NorthWest');
%
%\caption{Comparison in the analytical and numericall computed derivatives, 
% using MATLA \lstinline{gradient} function and the derivative matrix $\M{D}$ computed using
% the \lstinline{DOPbox}.}
%%
% Note there are significant error in the numerical compuation of the
% derivatives using the gradient function, i.e. with a 3 point approximation. 
% The result is grossly incorrect at the ends of the support. This
% indicates that the classical 3 term differential estimate is not suitable
% for use with boundary value problems.
% Zooming in on the plot, shows that there are significant errors over the 
% complete support in the derivatives computed using the 3 term estimate. 
%
% In contrast the derivatives computed using the matrix determined from the
% basis functions and their derivatives shows no significant errors. The
% matrix $\M{D}$ is well suited for the estimation of derivatives.
%
figure(fig5);
axis([-0.95,0.95,-0.02,0.01]);
% \caption{Errors in the estimated derivatives in the centre region of the support.}
%%
%
% The strength of the regularizing differentiating matrix is best seen when
% estimating derivatives in the presence of noise.
%
noiseGain = 0.05;
noise = noiseGain * norm(x) * randn( size(x));
yn = y + noise;
%
dyNMat = DyMat * yn;
dyNDop = Dr * yn;
%
fig5 = figure;
plot( x, dy - dyNMat, 'k');
hold on;
plot( x, dy - dyNDop, 'r');
grid on;
xlabel('Support');
ylabel('Error in $$d y(x)$$');
legend('Matlab Gradient','dop derivative','location','NorthWest');
% \caption{Estimation of derivatives in the presence of noise.}\label{fig:yn}
%%
% \subsection{Generate a matrix $\M{D}_L$ to perform local derivative approximation -- \lstinline{dopDiffLocal.m}}
%
% This function generates a local differentiating matrix. It supports both
% sparse and full matrix formats. The full matrix format is the default.
% Selecting the sparse option will generate a sparse matrix, this can be
% very important when dealing with large sparse sets of equations.
% 
% If we are looking to solve BVPs and IBVPs it will be necessary
% to solve problems of the form,
%
% \begin{equation}
%   \dot{\V{y}} = \M{D} \, \V{y}.
% \end{equation}
%
% Given measurements of $\dot{\V{y}}$ we may wist to compute $\V{y}$,
%
% \begin{equation}
%   \V{y} = \M{D}^+ \, \dot{\V{y}} + \nullS{\M{D}} \, \V{\alpha}.
% \end{equation}
%
% It is known that a differential operator as a single null, i.e. the
% constant vector $\V{1} \, \alpha$. Consequently, the matrix $\M{D}$ 
% should also have a single constant vector as a null space. This is not
% the case for the regularizing differential operator $\M{D}_r$. 
%
% Differentiating matrices generally become degenerate if computed globally 
% for a large number of points and a high degree, this is a fundamental 
% property of such matrices.
% Comsequently, the \lstinline{DOPbox} provides the function
% \lstinline{dopDiffLocal.m} to compute a local polynomial differential
% approximator.
%
ls = 7;
noBfs = ls;
%
Dl = dopDiffLocal(x, ls, noBfs);
dyL = Dl * y;
%
%
fig6 = figure;
plot( x, dy - dyDop, 'r');
hold on;
plot( x, dy - dyL, 'k');
grid on;
xlabel('Support');
ylabel('Error in $$d y(x)$$');
legend('dop derivative','local derivative','location','NorthWest');
%
%\caption{Comparison of the regularizing and local estimates for the derivative of $y(x)$}
%%
% There are indeed larger errors in the local derivative approximation, but
% at a scale of $10^{-12}$ there are very small. The matrix has, despite
% being dimension $100 \times 100$ the correct null space. 
%
% A numerical test for the quality of the matrix is the use of its pseudo
% inverse as a numerical integrator.
%
yMat = pinv(DyMat) * dy;
yDopLocal = pinv(Dl) * dy;
%
% ensure that the curve is mean free, this eliminates the need to deal with
% the constant of integration.
%
yt = y - mean(y);
fig7 = figure;
plot( x, y - yMat, 'k');
hold on;
plot( x, y - yDopLocal, 'r');
grid on;
legend('Matlab Gradient','dopLocal','location','NorthWest');
%
% \caption{Error in the curve reconstruction from its derivatives using 
% the pseudo inverse of the gradient operator in MATLAB and the 
% differentiating matrix computed using dopLocal.}
%%
% This result clearly shows that the differentiating matrix computed using
% \lstinline{dopDiffLocal.m} is well suited for inverse problems.
%%
% \subsection{Apply a set of constraints to a set of basis functions -- \lstinline{dopConstrain.m}}
%
% Boundary value problems, initial value problems and problems with inner
% constraints are characterized by a differential equation and a set of
% constraints which must be fulfilled by the solutions. The function
% \lstinline{dopConstrain} enables the application of constraints
% to a set of basis functions.
%
% \begin{figure}[H]
%    \centering
%    \includegraphics[width=13.5099cm]{Figures/cantileverFig.eps}
%    \caption{A simple cantilever with its associated constraints.}
% \end{figure}
%
% A simple cantilever provides a good example of a boundary value problem.
% The bending of the beam is described by the Euler-Bernouli equation and
% the following constraints need to be fulfilled, $y(0) = 0$, $\dot{y}(0) = 1$
% and $\ddot{y}(1) = 0$. We now demonstrate generating a set of addmissible
% functions for this problem using \lstinline{dopConstrain}. This is an
% interesting simple example since it requires both value and derivative
% constraints.
%
n = 100;
noBfs = 7;
%
xs = linspace(0,1,n )';
%
% And generate a set of basis functions
%
B = dop( xs, noBfs );
%%
% Now it is necessary to define the constraints. 
c1 = zeros( size(xs) );
%
% the first location in c corresponds to x = 0; settin a 1 at this location
% required c1' * y = 0, i.e. it implements the constraint
%
c1(1) = 1;
%
% the second constraint is a derivative constraint. First we must compute
% the derivative matrix and extract the first row which computes the
% derivative at the first point
%
ls = 5;
D = dopDiffLocal( xs, ls, ls );
d1 = D(1,:);
c2 = d1';
% the third constraint in on the end of the beam and is on the second
% deravitive.
%
D2 = D * D;
d21 = D2(end,:);
c3 = d21';
%
% Now conatinating the constraints
%
C = [c1, c2, c3];
%
% And applying the constraints to the basis functions B
% 
Bc = dopConstrain( C, B );
%%
% Having generated the basis functions we can now view them
%
noBcs = noBfs - rank(C);
figure;
plot( x, Bc(:,1), 'k');
hold on;
for k=2:(noBcs)
    plot( x, Bc(:,k), 'k');
end;
xlabel('Support');
ylabel('Amplitude');
grid on;
%
%\caption{Constrained basis functions for the cantilever example.}
%%
% The generated basis functions are also orthonormal. This makes them
% optimal for the propagation of Gaussian noise.
%
Gc = Bc' * Bc  
GcPerp = eye( noBcs ) - Gc
%%
% \subsection{Global and local polynomial approximation -- \lstinline{dopApproxLocal.m}}
%%
% In this section both global and local polynomial approximation is
% demonstrated. Global approximation referrs to applying the basis functions
% to the full length of the support. Given a set of noise data $\V{y}_n$
% located at the points $\V{x}$, global approximation is performed by
% generating a set of basis functions $\M{B}_x$ of the desired degree $d$ at the
% nodes defined by $\V{x}$. The approximation is performed by computing the
% projection onto the basis functions,
%
% \begin{equation}
%   \V{y}_{ag} = \M{B}_x \, \MT{B}_x \, \V{y}_n.
% \end{equation}
%
% This can also be looked at in a simular manner to  \lstinline{polyfit}
% and \lstinline{polyval}. The coefficients of the Gram polynomial
% corresponding to $\V{y}_n$ are computed as,
%
% \begin{equation}
%       \V{g} = \MT{B}_x \, \V{y}_n,
% \end{equation}
%
% this corresponds to \lstinline{polyfit}. The reconstruction of the signal is
% performed as,
%
% \begin{equation}
%   \V{y}_{ag} = \M{B}_x \, \V{g}.
% \end{equation}
% 
% corresponding to \lstinline{polyval}.
%
% There is clearly a direct relationship between the coefficients of the
% Gram polynomials and the geometric polynomials: this is obtained by
% looking at the reconstruction equation,
%
% \begin{equation}
%   \V{y}_{ag} = \M{B}_x \, \V{g} = \M{V} \, \V{p}
% \end{equation}
%
% whereby, $\M{V}$ is the Vandermonde matrix and $\V{p}$ would be the
% coefficients delivered by \lstinline{polyfit}. Consequently, the
% coefficinets $\V{p}$ can be computed from $\V{g}$ with the following
% MATLAB code,
%
% \begin{lstlisting}
%   p = V\(Bx * g). 
% \end{lstlisting}
%
% At this point it is important to note that in some applications it will
% be possible to compute a perfectly good approximation using the Gram
% polynomials; however, if the Vandermonde matrix $\M{V}$ has become
% degenerate it will not be possible to compute the corresponding geometric
% coefficients $\V{p}$.
%
% Using the data from Figure~\ref{fig:yn} to demonstrate fitting
degree = 6;
noBfs = degree + 1;
%
% Generat the basis functions
%
Bx = dop( x, noBfs );
%
% Compute the spectrum
%
g = Bx' * yn
%
% Reconstruct, i.e., compute the approximation. 
%
yag = Bx * g;
%
figure;
subplot(2,1,1);
plot( x, yn, 'k.');
hold on;
plot( x, yag, 'k');
grid on;
ylabel('Fit');
legend('$$y_n$$', '$$y_{ag}$$','Location','NorthEast');
%
subplot(2,1,2);
plot( x, yn - yag, 'k.');
grid on;
ylabel('$$\epsilon = y_n - y_{ag}$$');
xlabel('Support');
%%
% However, there are many examples where a global approximation will suffer
% from the Runge phenomena. In such cases local polynomial approximations
% may provide better approximations. Consider the data shown in
% Figure~\ref{fig:gl} this is a data set for which Runge phenomena is
% relevant
%
load approxLocalData
%
figGL = figure;
plot( x, y, 'k');
hold on;
plot( x, yn, 'k.');
xlabel('Support');
ylabel('Value');
axis([-1,1,-0.1,1.1]);
grid on;
%
% \caption{Test data set to compare global and local polynomial
% approximation.}\label{fig:gl}
%%
% Approximating the above data set with a degree $d = 25$ polynomial
%
degree = 25;
B = dop( x, degree + 1 );
yg = B * B' * yn;
%
figure(figGL);
plot( x, yg, 'r');
legend('$$y(x)$$','$$y_n$$','$$y_g$$');
%
%\caption{The degree $d=25$ polynomial approximation exhibits a significant 
% effect from the Runga phenomena.}
%%
% 
% Now comparing the results to a local polynomial approximation (LPA). The 
% the LPA is characterized by the length of the support $l_s$, i.e. the
% number of points used for the local approximation and the degree $d$ of
% the approximation. The function \lstinline{dopApproxLocal.m}
% generates a matrix $\M{S}$ which performs the LPA.
%
ls = 11;
d = 3;
%
% Generate the linear operator
%
S = dopApproxLocal( x, ls, d + 1 );
%
% Compute the approximation
%
ys = S * y;
%
figGL2 = figure;
plot( x, y, 'k');
hold on;
plot( x, yn, 'k.');
xlabel('Support');
ylabel('Value');
axis([-1,1,-0.1,1.1]);
grid on;
plot( x, ys, 'r');
legend('$$y(x)$$','$$y_n$$','$$y_s$$');
%
% \caption{local polynomial approximation.}
%
%%
% \subsection{Interpolating using basis functions -- \lstinline{dopInterpolate.m}}
%
% The methods provided here generate basis functions for arbitrary nodes
% within the unit circle on the complex plane. In many cases the corresponding 
% analytical functions are not directly available. The function \lstinline{dopInterpolate}
% provided a measns of interpolating within such functions.
%
% The test data set \lstinline{interpolateData} cantains a set of
% non-uniformly spaced x values (xnu), a uniformly spaced set of nodes (x) 
% and the function y(x) (y) evaluated at these points. the
% example here is to reevaluate the equations are regularly spaced nodes.
%
load interpolateData;
%
% Generate a set of basis functions for original set of nodes
%
noBfs = 5;
[B, dB, rCfs] = dop( x, noBfs );
%
% compute the spectrum of y with respect to B.
%
polyCoeffs = B' * y;
yb = B * polyCoeffs;
%
% Interpolate onto the non-uniform set of nodes
%
yi = dopInterpolate( polyCoeffs, rCfs, xnu );
%
figure;
plot( x, y, 'k');
hold on;
plot( xnu, yi, 'ro','markerFaceColor', 'w');
plot( x, y, 'k.');
grid on;
xlabel('Support');
ylabel('Amplitude');
%\caption{Example of using \lstinline{dopInterpolate} to recompute a function at a new set of arbitrary nodes.}
%%
% \subsection{Generate a set sine basis functions -- \lstinline{dopSine.m}}
%
% Generate a set of sine wave basis functions and their derivatives.
%
nrS = 300;
nrBfs = 4;
[S, dS, x] = dopSine( nrS, nrBfs );
%
figS = figure;
subplot(2,1,1);
for k=1:nrBfs
    plot( x, S(:,k), 'k');
    hold on;
end;
grid on;
ylabel('sine(x)');
subplot(2,1,2);
for k=1:nrBfs
    plot( x, dS(:,k), 'k');
    hold on;
end;
grid on;
ylabel('d sine(x) / dx');
xlabel('x');
%
% \caption{Sine basis functions and their derivatives.}
%%
% \subsection{Discrete Weighted Orthogonal Polynomials: \lstinline{dopWeightBasis}}
%
% Discrete orthogonal weighted polynomials are defined on the orthogonality condition,
% \begin{equation}
%   \M{B} \, \M{W} \, \MT{B} = \M{I}
% \end{equation}
%
% Whereby, $\M{W}$ is a full rank positive definite weighting matrix,
% it is not restricted to being a diagonal matrix. 

%
% Example generate 10 weighted basis functions
%
nw = 100;
nrBfs = 10;
%
% generate the Gram nodes
%
xw = dopNodes( nw, 'gram');
%
% generate the Unitary basis functions
%
B = dop( xw , nrBfs );
%
% Define the weights and compute the weighted basis functions. These
% correspond to the Gegenbauer polynomials for alpha = -0.5 and scaled to
% be unitary in the definition B' W B = I.
%
w = 1./(1 - xw.^2);
W = diag( w );
Bw = dopWeightBasis( B, W ); 
%
figW = figure;
plot( xw, Bw(:,2), 'k');
hold on;
plot( xw, Bw(:,3), 'k');
plot( xw, Bw(:,4), 'k');
plot( xw, B(:,2), 'r');
plot( xw, B(:,3), 'r');
plot( xw, B(:,4), 'r');
xlabel('Support');
ylabel('Weighted polynomial value');
grid on;
%
% \caption{Discrete othogonal weighted polynomials $\MT{B}\,\M{W}\,\M{B} = \M{I}$. 
% The wighting used here is $w = \frac{1}{1 - x^2}$. Weighted polynomials in black and unweighted in red.}
%%
% \subsection{Weighted regression with weighted basis functions}
%
% The weighted basis functions presented in the previous section can used
% for a simple implementation of a weighted regression, e.g., a covariance
% weighted regression.
%
% The cost function for weigthed regression is,
% %
% \begin{equation}
%   E = \left(\V{y} - \M{B}_w \, \V{\alpha}^\mathrm{T} \right) \, \M{W}
%   \left(\V{y} - \M{B}_w \, \V{\alpha} \right).
% \end{equation}
% %
% In the case of covariance weighted regression $\M{W} = \V{\Lambda^{-1}}$
% whereby $\Lambda^{-1}$ is the symmetric positive definite covariance matrix.
% Now expanding the above equation and differentiating with respect
% to $\V{\alpha}$ yields,
% %
% \begin{equation}
%   \frac{\md E}{\md \V{\alpha}} =
%       - 2 \, \MT{B}_w \, \M{W} \, \V{y} + 2 \, \MT{B}_w\, \M{W} \,\M{B}_w 
%       \, \V{\alpha}.
% \end{equation}
% %
% Evaluating the derivative equal to the zero vector, yields the so-called
% normal equations for weighted regression,
% %
% \begin{equation}
%   \MT{B}_w\, \M{W} \,\M{B}_w \, \V{\alpha} = \MT{B}_w \, \M{W} \, \V{y}.
% \end{equation}
% %
% If the weighted basis functions are generated such that $\MT{B}_w\, \M{W}
% \,\M{B}_w = \M{I}$ as shown in the previous section. 
% Then the normal equations simplify to:
% %
% \begin{equation}
%   \V{\alpha} = \MT{B}_w \, \M{W} \, \V{y}.
% \end{equation}
% %
% Computing in this manner avoide the necessity to calculate matrix inverses
% at the time of fitting. This greatly simplifies the computation 
% of weighted regression.
%
%%
% \subsection{\lstinline{dopFit} and \lstinline{dopVal}}
%
% These are two wrapper functions to simplify the use of the DOP library
% for fitting and evaluating discrete orthogonal polynomials. They are the
% discrete orthogonal polynomial equivalent of \lstinline{polyfit} and 
% \lstinline{polyval}. However, the discrete orthogonal polynomials can
% perform fitting up to degree of $D = 1000$ without serious numerical
% errors. 
%
% The function \lstinline{dopVit} also returen the covariance matrix for
% the polynomial coefficients and \lstinline{dopVal} returen the covariance 
% matrix for the fitted $\V{y}_f$ values. 
%
% A Taylor expansion for a function is equivalent to computing a polynomial
% approximation. The discrete orthogonal polynomials enable the computation
% of a Gram expansion, i.e. a polynomial expansion, of very high degree. 
% This may be advantegous, e.g., when solving differential equations which
% are known to have no analytical solutions such as the Mathieu
% differential equation.
%
%%
% \subsubsection{Fitting a high degree polynomials}
%
% In this first test a polynomial of degree $d=100$ is fit to a data set.
% Fitting with this high degree is not possible with the standard MATLAB
% function \lstinline{polyfit}.
%
% Load a test data set:
load dopFitData;
%
% Perform the fit and evaluate the polynomial
%
d = 100;
[g, Lg, S] = dopFit( x, y, d );
yg = dopVal( g, S );
%
% Compute the residual
%
r = y - yg;
%
% Plot the fit
%
fitFig1 = figure;
plot( x, y, 'k.');
hold on;
plot( x, yg, 'r');
ylabel('Value');
grid on;
xlabel('x');
%
%\caption{Example of a discrete orthogonal polynomial fit with degree $d=100$.}
%
%%
% Plot the residual
%
fitFig2 = figure;
plot( x, r, 'k.');
ylabel('Residual');
grid on;
%
%\caption{Residual of the discrete orthogonal polynomial fit with degree $d=100$.}
%
%%
% Compute the histogram of the residual
%
nrBins = 15;
fitFig3 = figure;
hist( r, nrBins );
xlabel('Residual');
ylabel('Frequency');
%
%\caption{Histogram of the residual of the discrete orthogonal polynomial 
% fit with degree $d=100$. This demonstrated the Gaussian nature of the 
% residual, i.e. the fit has not produced systematic errors.}
%
%%
% \subsubsection{High degree polynomial expansion}
%
% This section demonstrates the application of \lstinline{dopFit} to
% compute a high degree polynomial expansion for a function. The cosine
% finction is chosen here because its Taylor expansion is known;
% furthermore, \lstinline{polyfit} fails to model this data correctly.
%
% Generate the test function
%
%
nrCycles = 15;
nrPts = 300;
%
x = dopNodes( nrPts, 'GramEnds');
%
phi = linspace( 0, 2*pi*nrCycles, nrPts )';
y = cos(phi);
%
% compute the fit
%
d = 70;
[g, Lg, S] = dopFit( x, y, d );
%
% Note the structure S also delivers the values of y fit. These value are
% required internally to do the covariance computations.
%
expFig = figure;
plot(x, y, 'k.');
hold on;
plot( x, S.yg, 'r');
xlabel('x');
ylabel('Value');
%
%\caption{A discrete orthogonal polynomial expansion of degree $d=70$ for a
% cosine. This demonstrated the possability of computing a polynomial
% expandions for functions which require a high degree.}
%
%%
% \bibliographystyle{plain}
% \bibliography{harkerOleary}