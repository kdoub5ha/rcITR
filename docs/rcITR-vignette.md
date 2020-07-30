<!DOCTYPE html>

<html xmlns="http://www.w3.org/1999/xhtml">

<head>

<meta charset="utf-8" />
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="pandoc" />

<meta name="viewport" content="width=device-width, initial-scale=1">



<title>rcITR-vignette</title>



<style type="text/css">code{white-space: pre;}</style>
<style type="text/css" data-origin="pandoc">
a.sourceLine { display: inline-block; line-height: 1.25; }
a.sourceLine { pointer-events: none; color: inherit; text-decoration: inherit; }
a.sourceLine:empty { height: 1.2em; }
.sourceCode { overflow: visible; }
code.sourceCode { white-space: pre; position: relative; }
div.sourceCode { margin: 1em 0; }
pre.sourceCode { margin: 0; }
@media screen {
div.sourceCode { overflow: auto; }
}
@media print {
code.sourceCode { white-space: pre-wrap; }
a.sourceLine { text-indent: -1em; padding-left: 1em; }
}
pre.numberSource a.sourceLine
  { position: relative; left: -4em; }
pre.numberSource a.sourceLine::before
  { content: attr(title);
    position: relative; left: -1em; text-align: right; vertical-align: baseline;
    border: none; pointer-events: all; display: inline-block;
    -webkit-touch-callout: none; -webkit-user-select: none;
    -khtml-user-select: none; -moz-user-select: none;
    -ms-user-select: none; user-select: none;
    padding: 0 4px; width: 4em;
    color: #aaaaaa;
  }
pre.numberSource { margin-left: 3em; border-left: 1px solid #aaaaaa;  padding-left: 4px; }
div.sourceCode
  {  }
@media screen {
a.sourceLine::before { text-decoration: underline; }
}
code span.al { color: #ff0000; font-weight: bold; } /* Alert */
code span.an { color: #60a0b0; font-weight: bold; font-style: italic; } /* Annotation */
code span.at { color: #7d9029; } /* Attribute */
code span.bn { color: #40a070; } /* BaseN */
code span.bu { } /* BuiltIn */
code span.cf { color: #007020; font-weight: bold; } /* ControlFlow */
code span.ch { color: #4070a0; } /* Char */
code span.cn { color: #880000; } /* Constant */
code span.co { color: #60a0b0; font-style: italic; } /* Comment */
code span.cv { color: #60a0b0; font-weight: bold; font-style: italic; } /* CommentVar */
code span.do { color: #ba2121; font-style: italic; } /* Documentation */
code span.dt { color: #902000; } /* DataType */
code span.dv { color: #40a070; } /* DecVal */
code span.er { color: #ff0000; font-weight: bold; } /* Error */
code span.ex { } /* Extension */
code span.fl { color: #40a070; } /* Float */
code span.fu { color: #06287e; } /* Function */
code span.im { } /* Import */
code span.in { color: #60a0b0; font-weight: bold; font-style: italic; } /* Information */
code span.kw { color: #007020; font-weight: bold; } /* Keyword */
code span.op { color: #666666; } /* Operator */
code span.ot { color: #007020; } /* Other */
code span.pp { color: #bc7a00; } /* Preprocessor */
code span.sc { color: #4070a0; } /* SpecialChar */
code span.ss { color: #bb6688; } /* SpecialString */
code span.st { color: #4070a0; } /* String */
code span.va { color: #19177c; } /* Variable */
code span.vs { color: #4070a0; } /* VerbatimString */
code span.wa { color: #60a0b0; font-weight: bold; font-style: italic; } /* Warning */

</style>
<script>
// apply pandoc div.sourceCode style to pre.sourceCode instead
(function() {
  var sheets = document.styleSheets;
  for (var i = 0; i < sheets.length; i++) {
    if (sheets[i].ownerNode.dataset["origin"] !== "pandoc") continue;
    try { var rules = sheets[i].cssRules; } catch (e) { continue; }
    for (var j = 0; j < rules.length; j++) {
      var rule = rules[j];
      // check if there is a div.sourceCode rule
      if (rule.type !== rule.STYLE_RULE || rule.selectorText !== "div.sourceCode") continue;
      var style = rule.style.cssText;
      // check if color or background-color is set
      if (rule.style.color === '' || rule.style.backgroundColor === '') continue;
      // replace div.sourceCode by a pre.sourceCode rule
      sheets[i].deleteRule(j);
      sheets[i].insertRule('pre.sourceCode{' + style + '}', j);
    }
  }
})();
</script>



<style type="text/css">body {
background-color: #fff;
margin: 1em auto;
max-width: 700px;
overflow: visible;
padding-left: 2em;
padding-right: 2em;
font-family: "Open Sans", "Helvetica Neue", Helvetica, Arial, sans-serif;
font-size: 14px;
line-height: 1.35;
}
#header {
text-align: center;
}
#TOC {
clear: both;
margin: 0 0 10px 10px;
padding: 4px;
width: 400px;
border: 1px solid #CCCCCC;
border-radius: 5px;
background-color: #f6f6f6;
font-size: 13px;
line-height: 1.3;
}
#TOC .toctitle {
font-weight: bold;
font-size: 15px;
margin-left: 5px;
}
#TOC ul {
padding-left: 40px;
margin-left: -1.5em;
margin-top: 5px;
margin-bottom: 5px;
}
#TOC ul ul {
margin-left: -2em;
}
#TOC li {
line-height: 16px;
}
table {
margin: 1em auto;
border-width: 1px;
border-color: #DDDDDD;
border-style: outset;
border-collapse: collapse;
}
table th {
border-width: 2px;
padding: 5px;
border-style: inset;
}
table td {
border-width: 1px;
border-style: inset;
line-height: 18px;
padding: 5px 5px;
}
table, table th, table td {
border-left-style: none;
border-right-style: none;
}
table thead, table tr.even {
background-color: #f7f7f7;
}
p {
margin: 0.5em 0;
}
blockquote {
background-color: #f6f6f6;
padding: 0.25em 0.75em;
}
hr {
border-style: solid;
border: none;
border-top: 1px solid #777;
margin: 28px 0;
}
dl {
margin-left: 0;
}
dl dd {
margin-bottom: 13px;
margin-left: 13px;
}
dl dt {
font-weight: bold;
}
ul {
margin-top: 0;
}
ul li {
list-style: circle outside;
}
ul ul {
margin-bottom: 0;
}
pre, code {
background-color: #f7f7f7;
border-radius: 3px;
color: #333;
white-space: pre-wrap; 
}
pre {
border-radius: 3px;
margin: 5px 0px 10px 0px;
padding: 10px;
}
pre:not([class]) {
background-color: #f7f7f7;
}
code {
font-family: Consolas, Monaco, 'Courier New', monospace;
font-size: 85%;
}
p > code, li > code {
padding: 2px 0px;
}
div.figure {
text-align: center;
}
img {
background-color: #FFFFFF;
padding: 2px;
border: 1px solid #DDDDDD;
border-radius: 3px;
border: 1px solid #CCCCCC;
margin: 0 5px;
}
h1 {
margin-top: 0;
font-size: 35px;
line-height: 40px;
}
h2 {
border-bottom: 4px solid #f7f7f7;
padding-top: 10px;
padding-bottom: 2px;
font-size: 145%;
}
h3 {
border-bottom: 2px solid #f7f7f7;
padding-top: 10px;
font-size: 120%;
}
h4 {
border-bottom: 1px solid #f7f7f7;
margin-left: 8px;
font-size: 105%;
}
h5, h6 {
border-bottom: 1px solid #ccc;
font-size: 105%;
}
a {
color: #0033dd;
text-decoration: none;
}
a:hover {
color: #6666ff; }
a:visited {
color: #800080; }
a:visited:hover {
color: #BB00BB; }
a[href^="http:"] {
text-decoration: underline; }
a[href^="https:"] {
text-decoration: underline; }

code > span.kw { color: #555; font-weight: bold; } 
code > span.dt { color: #902000; } 
code > span.dv { color: #40a070; } 
code > span.bn { color: #d14; } 
code > span.fl { color: #d14; } 
code > span.ch { color: #d14; } 
code > span.st { color: #d14; } 
code > span.co { color: #888888; font-style: italic; } 
code > span.ot { color: #007020; } 
code > span.al { color: #ff0000; font-weight: bold; } 
code > span.fu { color: #900; font-weight: bold; }  code > span.er { color: #a61717; background-color: #e3d2d2; } 
</style>




</head>

<body>




<h1 class="title toc-ignore">rcITR-vignette</h1>



<div id="introduction" class="section level1">
<h1>Introduction</h1>
<p>The package <code>rcITR</code> is an individualized treatment rule (ITR) discovery tool that estimates ITRs under a constrained optimization framework. The main functions are <code>rcDT()</code> and <code>rcRF()</code>. The <code>rcDT()</code> function constructs a risk controlled decision tree (rcDT) and <code>rcRF()</code> constructs a risk controlled random forest (rcRF). The aim of each model is to maximize an expected benefit (Y) while constraining risk (R) at a pre-determined threshold. Each model is constructed to optimize a Lagrangian objective function which takes the form, <span class="math display">\[
L(d) = E(Y) - \lambda(E(R) - \tau).
\]</span> Here, <span class="math inline">\(d\)</span> is a treatment decision rule of the form <span class="math inline">\(d(X): X \rightarrow A\)</span> for predictors <span class="math inline">\(X\)</span> and treatments <span class="math inline">\(A\in \{0,1\}\)</span>. <span class="math inline">\(E(Y)\)</span> and <span class="math inline">\(E(R)\)</span> are estimated using the so-called “value” function and correspond to the expected benefit and risk if rule <span class="math inline">\(d\)</span> were used. For example, the estimated value of treatment recommendation <span class="math inline">\(d\)</span> would be <span class="math display">\[E(Y) =N^{-1} \sum_{i=1}^{N} \frac{\textbf{I}_{a_i=d(\textbf{x}_i)}}{\hat{p}(a_i|\textbf{x}_i)}y_i \]</span> where <span class="math inline">\(\textbf{x}_i\)</span> is the covariate vector, <span class="math inline">\(a_i\)</span> is the original treatment assignment, <span class="math inline">\(d(\textbf{x}_i)\)</span> is the proposed treatment rule, and <span class="math inline">\(\hat{p}(a_i|\textbf{x}_i)\)</span> is the probability of receiving treatment given the covariates. The estimated ``value&quot; associated with the risk outcome is simiarly estimated. The pre-determined risk threshold <span class="math inline">\(\tau\)</span> is set with knowledge of the expected range of risk values in the population of interest and <span class="math inline">\(\lambda\)</span> is the penalty for rule <span class="math inline">\(d\)</span> exceeding the risk threshold.</p>
<p>The input dataset needs to be a data.frame() with the following: (1) benefit score column, preferably labeled as <code>y</code>, (2) risk score column, preferably labeled as <code>r</code>, (3) binary treatment column, preferably labeled <code>trt</code> having treatment indicator scores of 0 or +1, (4) propensity score column, preferably labeled <code>prtx</code>, (5) columns of predictors</p>
</div>
<div id="constructing-an-rcdt-model" class="section level1">
<h1>Constructing an rcDT Model</h1>
<p>A single rcDT tree structrue can be constructed using the function <code>rcDT</code>. By default the number of observations allowed in a terminal node is 20 (<code>min.ndsz = 20</code>), there must be at least 5 observations from each treatment group in a terminal node (<code>n0 = 5</code>), and the maximum tree depth is set at 15 (<code>max.depth = 15</code>). The initial objective value of the root node is the maximum of <span class="math inline">\(\hat{L}(d)\)</span> with all subjects given treatment or all patients given control. The root node is split only made if there is a binary partition on a predictor for which the Lagrangian objective increases. The same is true for additional splits so that the tree cannot grow larger unless there is an increase in overall objective value of the tree given by the split. The function requires the user to input values for <span class="math inline">\(\tau\)</span> and <span class="math inline">\(\lambda\)</span>.</p>

<div id="rcdt-example" class="section level2">
<h2>rcDT Example</h2>
<p>We will use an example dataset generated from the function <code>generateData()</code> which simulates RCT data from the model with</p>
<p><span class="math display">\[ Y = 1 - X_3 + X_4 + 3\text{I}_{x \in S}\text{I}_{A = 1} + 3\text{I}_{x \notin S}\text{I}_{A = 0} + \epsilon_Y. \]</span> <span class="math display">\[ R = 2 + 2X_3 + X_4 + (2X_2 + 0.1)(2A-1) + \epsilon_R. \]</span> where <span class="math inline">\(S = \{X_1\le 0.7\hspace{2mm} \cap\hspace{2mm} X_2\le 0.7\}\)</span>, <span class="math inline">\(\epsilon_Y\sim N(0,1)\)</span>, and <span class="math inline">\(\epsilon_R\sim N(0,0.5)\)</span>. Covariates <span class="math inline">\(X_1 - X_{10} \sim \text{Unif}(0,1)\)</span>. The simulated data has 1000 observations.</p>

<div class="sourceCode" id="cb1"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb1-1" title="1"><span class="kw">set.seed</span>(<span class="dv">123</span>)</a>
<a class="sourceLine" id="cb1-2" title="2">dat &lt;-<span class="st"> </span><span class="kw">generateData</span>(<span class="dt">n =</span> <span class="dv">1000</span>)</a>
<a class="sourceLine" id="cb1-3" title="3"><span class="kw">str</span>(dat)</a></code></pre></div>
<pre><code>## &#39;data.frame&#39;:    1000 obs. of  14 variables:
##  $ x1  : num  0.274 0.594 0.16 0.853 0.848 0.478 0.774 0.295 0.066 0.441 ...
##  $ x2  : num  0.16 0.145 0.149 0.514 0.493 0.616 0.447 0.056 0.005 0.222 ...
##  $ x3  : num  0.206 0.943 0.379 0.626 0.184 0.659 0.735 0.496 0.869 0.832 ...
##  $ x4  : num  0.304 0.833 0.594 0.807 0.294 0.141 0.889 0.008 0.569 0.968 ...
##  $ x5  : num  0.249 0.989 0.717 0.652 0.241 0.087 0.297 0.542 0.781 0.272 ...
##  $ x6  : num  0.44 0.397 0.372 0.529 0.074 0.717 0.243 0.844 0.995 0.105 ...
##  $ x7  : num  0.93 0.591 0.08 0.845 0.03 0.815 0.391 0.143 0.667 0.619 ...
##  $ x8  : num  0.578 0.489 0.742 0.226 0.749 0.825 0.1 0.276 0.021 0.596 ...
##  $ x9  : num  0.798 0.879 0.242 0.56 0.905 0.008 0.163 0.26 0.773 0.071 ...
##  $ x10 : num  0.311 0.325 0.87 0.329 0.126 0.356 0.931 0.875 0.82 0.022 ...
##  $ trt : num  0 1 0 1 1 0 1 1 1 0 ...
##  $ y   : num  -0.259 2.597 -0.302 2.04 -0.105 ...
##  $ r   : num  2.44 5.06 3.4 5.58 3.56 ...
##  $ prtx: num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...</code></pre>

<p>Figure 1 presents the risk and efficacy score distributions for the original treated and control groups. Risk scores are lower on average for the control arm of the study and so we may want to set the risk threshold to be around the center of the distribution of risk scores for the control arm. Note that efficacy scores do not differ much between the two groups. We will set <span class="math inline">\(\tau = 2.75\)</span>, and use a penalty of <span class="math inline">\(\lambda = 1\)</span> for the demonstration.</p>
<div class="figure">
<img src="data:image/png;base64,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" alt="Figure 1. Risk score distribution for simulated data" />
<p class="caption">Figure 1. Risk score distribution for simulated data</p>
</div>

<p>Below is the tree structure generated for the simulated data corresponding the a risk threshold of <span class="math inline">\(\tau = 2.75\)</span> and penalty <span class="math inline">\(\lambda = 1\)</span>.</p>
<div class="sourceCode" id="cb3"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb3-1" title="1">tre1 &lt;-<span class="st"> </span><span class="kw">rcDT</span>(<span class="dt">data =</span> dat, </a>
<a class="sourceLine" id="cb3-2" title="2">             <span class="dt">split.var =</span> <span class="dv">1</span><span class="op">:</span><span class="dv">10</span>,</a>
<a class="sourceLine" id="cb3-3" title="3">             <span class="dt">risk.threshold =</span> <span class="fl">2.75</span>,</a>
<a class="sourceLine" id="cb3-4" title="4">             <span class="dt">lambda =</span> <span class="dv">1</span>,</a>
<a class="sourceLine" id="cb3-5" title="5">             <span class="dt">efficacy =</span> <span class="st">&quot;y&quot;</span>,</a>
<a class="sourceLine" id="cb3-6" title="6">             <span class="dt">risk =</span> <span class="st">&quot;r&quot;</span>,</a>
<a class="sourceLine" id="cb3-7" title="7">             <span class="dt">col.trt =</span> <span class="st">&quot;trt&quot;</span>,</a>
<a class="sourceLine" id="cb3-8" title="8">             <span class="dt">col.prtx =</span> <span class="st">&quot;prtx&quot;</span>)</a>
<a class="sourceLine" id="cb3-9" title="9">tre1<span class="op">$</span>tree</a></code></pre></div>
<pre><code>##   node size n.1 n.0  trt.effect var vname cut.1 cut.2     score
## 1    0 1000 493 507 -0.06270569   2    x2     l 0.288 0.6925219
## 2   01  304 139 165  1.65841703   1    x1     l   0.7 0.9232185
## 4  011  223 104 119  3.16167696  NA  &lt;NA&gt;  &lt;NA&gt;  &lt;NA&gt;        NA
## 5  012   81  35  46 -2.52519851  NA  &lt;NA&gt;  &lt;NA&gt;  &lt;NA&gt;        NA
## 3   02  696 354 342 -0.81243281  NA  &lt;NA&gt;  &lt;NA&gt;  &lt;NA&gt;        NA</code></pre>
<p>The output contains a summary of the tree structure. The <code>node</code> column begins with the root node <code>0</code> and each subsequent number indicates the direction of the split, with <code>1</code> indicating the left (less than or equal to) node and <code>2</code> indicating the right (greater than) node. The first row indicates that the covariate <span class="math inline">\(X_2\)</span> is selected as the first splitting variable with a cut point of <code>cut.2 = 0.375</code>. The decision is to send treatment to the left node (<code>cut.1 = &quot;l&quot;</code>). <code>size</code>, <code>n.1</code>, and <code>n.0</code> indicate there are observations in the root node, with originally treated and originally on control. The second row with <code>node = 01</code> contains information from the left child node with interpretations similar to those described for the root node. The splitting information denoted by <code>NA</code> indicates a terminal node, <code>01111</code> for instance.</p>

<p>If we use the same risk threshold of <span class="math inline">\(\tau = 2.75\)</span> but with a stricter penalty <span class="math inline">\(\lambda = 2\)</span>, we obtain the following tree structre.</p>
<div class="sourceCode" id="cb5"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb5-1" title="1">tre2 &lt;-<span class="st"> </span><span class="kw">rcDT</span>(<span class="dt">data =</span> dat, </a>
<a class="sourceLine" id="cb5-2" title="2">             <span class="dt">split.var =</span> <span class="dv">1</span><span class="op">:</span><span class="dv">10</span>,</a>
<a class="sourceLine" id="cb5-3" title="3">             <span class="dt">risk.threshold =</span> <span class="fl">2.75</span>,</a>
<a class="sourceLine" id="cb5-4" title="4">             <span class="dt">lambda =</span> <span class="dv">2</span>,</a>
<a class="sourceLine" id="cb5-5" title="5">             <span class="dt">efficacy =</span> <span class="st">&quot;y&quot;</span>,</a>
<a class="sourceLine" id="cb5-6" title="6">             <span class="dt">risk =</span> <span class="st">&quot;r&quot;</span>,</a>
<a class="sourceLine" id="cb5-7" title="7">             <span class="dt">col.trt =</span> <span class="st">&quot;trt&quot;</span>,</a>
<a class="sourceLine" id="cb5-8" title="8">             <span class="dt">col.prtx =</span> <span class="st">&quot;prtx&quot;</span>)</a>
<a class="sourceLine" id="cb5-9" title="9">tre2<span class="op">$</span>tree</a></code></pre></div>
<pre><code>##   node size n.1 n.0  trt.effect var vname cut.1 cut.2     score
## 1    0 1000 493 507 -0.06270569   2    x2     l 0.222 0.8998027
## 2   01  239 105 134  1.64864910   1    x1     l   0.7 1.0689384
## 5  011  178  81  97  3.17618787  NA  &lt;NA&gt;  &lt;NA&gt;  &lt;NA&gt;        NA
## 4  012   61  24  37 -2.94499308  NA  &lt;NA&gt;  &lt;NA&gt;  &lt;NA&gt;        NA
## 3   02  761 388 373 -0.60015177  NA  &lt;NA&gt;  &lt;NA&gt;  &lt;NA&gt;        NA</code></pre>
<p>Note that the tree structures generated for different penalty values are very different. The estimated risk in the training set when <span class="math inline">\(\lambda = 1\)</span> is 2.53 and when <span class="math inline">\(\lambda = 2\)</span> the estimated risk is 2.47. Thus we see that selection of <span class="math inline">\(\lambda\)</span> is important as both penalties controls risk below <span class="math inline">\(\tau = 2.75\)</span>, but setting <span class="math inline">\(\lambda = 2\)</span> may over control risk and result a loss of potential benefit. </p>
</div>
<div id="pruning-a-tree" class="section level2">
<h2>Pruning a Tree</h2>
<p>The function <code>prune()</code> determines the series of optimally pruned subtrees from a large tree. A branch of the large tree is deemed the “weakest” branch if trimming the branch from the large tree results in the smallest decrease in objective function value. This is repeated until the null tree is reached. Below is the result from the pruning procedure for <span class="math inline">\(\tau = 2.75\)</span> and <span class="math inline">\(\lambda = 1\)</span>. Each row corresponds to the sequence of pruning to be performed in the input tree. For instance, in this case we would first trim off descendant nodes of 0111 and declare <code>0111</code> a terminal node. The column <code>V</code> gives the sequence of objective function values for each pruning event.</p>

<div class="sourceCode" id="cb7"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb7-1" title="1">pruned1 &lt;-<span class="st"> </span><span class="kw">prune</span>(<span class="dt">tre =</span> tre1, </a>
<a class="sourceLine" id="cb7-2" title="2">                 <span class="dt">a =</span> <span class="dv">0</span>, </a>
<a class="sourceLine" id="cb7-3" title="3">                 <span class="dt">risk.threshold =</span> <span class="fl">2.75</span>, </a>
<a class="sourceLine" id="cb7-4" title="4">                 <span class="dt">lambda =</span> <span class="dv">1</span>)</a></code></pre></div>
<pre><code>##   subtree node.rm size.tree size.tmnl    alpha     V Benefit  Risk
## 1       1      01         5         3    0.231 0.923   0.747 2.574
## 2       2       0         3         2    0.383 0.693   0.555 2.612
## 3       3      NA         1         1 9999.000 0.310   0.031 2.471</code></pre>

<p>The first row represents the entire tree (<code>subtree 1</code>) and provides the weakest link (<code>node.rm</code>), number of total nodes in the subtree (<code>size.tree</code>), number of terminal nodes in the subtree (<code>size.tmnl</code>), the <span class="math inline">\(\alpha\)</span> parameter, the Lagrangian “value” for the subtree (<code>V</code>), the benefit in the training data (<code>Benefit</code>), and the risk (<code>Risk</code>). Here, the risk is controlled in all subtrees and the third subtree also gives the greatest benefit.</p>
</div>
<div id="cross-validation-for-model-selection" class="section level2">
<h2>Cross Validation for Model Selection</h2>
<p>One issue with the approach above for model selection is the risk of overfitting through using the training data alone for model selection. The function <code>rcDT.select()</code> will perform n-fold cross validation to select the optimal tuning parameter (<span class="math inline">\(\alpha\)</span>). The function returns the optimal model, selected tuning parameter, and several summary measures.</p>
<div class="sourceCode" id="cb9"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb9-1" title="1">rcDT.fit &lt;-<span class="st"> </span><span class="kw">rcDT.select</span>(<span class="dt">data =</span> dat, </a>
<a class="sourceLine" id="cb9-2" title="2">                        <span class="dt">split.var =</span> <span class="dv">1</span><span class="op">:</span><span class="dv">10</span>,</a>
<a class="sourceLine" id="cb9-3" title="3">                        <span class="dt">risk.threshold =</span> <span class="fl">2.75</span>, </a>
<a class="sourceLine" id="cb9-4" title="4">                        <span class="dt">efficacy =</span> <span class="st">&quot;y&quot;</span>,</a>
<a class="sourceLine" id="cb9-5" title="5">                        <span class="dt">risk =</span> <span class="st">&quot;r&quot;</span>,</a>
<a class="sourceLine" id="cb9-6" title="6">                        <span class="dt">col.trt =</span> <span class="st">&quot;trt&quot;</span>,</a>
<a class="sourceLine" id="cb9-7" title="7">                        <span class="dt">col.prtx =</span> <span class="st">&quot;prtx&quot;</span>,</a>
<a class="sourceLine" id="cb9-8" title="8">                        <span class="dt">nfolds =</span> <span class="dv">5</span>,</a>
<a class="sourceLine" id="cb9-9" title="9">                        <span class="dt">verbose =</span> <span class="ot">FALSE</span>)</a></code></pre></div>
<p>Note that the optimal tree presented here is selected using 5-fold cross validation and is <code>subtree 1</code> from the pruning procedure. The final model has a training set benefit estimated as 3.53 and risk estimated as 2.71.</p>
</div>
<div id="constructing-a-risk-controlled-random-forest-rcrf-model" class="section level2">
<h2>Constructing a Risk Controlled Random Forest (rcRF) Model</h2>
<p>A single tree which is trained using all available data may be overfitted, having poor predictive ability in an external dataset. Hence, we make a decision rule using a forest of rcDT predictors. Each rcDT model is constructed to be more variable, but the aggregation of the trees in the mitigates this variance. An rcRF is contructed using the function <code>rcRF()</code> and requires similar inputs to the <code>rcDT()</code> function. To randomized the growth of trees in the forest a subset of predictors, <code>mtry</code>, is selected as potential splitting variables at each split which defaults to the maximum of 1/3 the number of splitting variables and 1. By default the number of trees contructed, <code>ntree</code>, is 500. Each tree is grown using a bootstrap sample taken from the input dataset without replacement. The function returns the bootstrap samples used in tree construction, the trees, the model parameters, and several summaries for the in- and out-of-bag samples. We leverage the out-of-bag sample, i.e. observations not in the bootstrap sample, to obtain risk control estimates which are less biased than those obtained directly from training data. Hence, several values of <span class="math inline">\(\lambda\)</span> may need to be used in order to final an optimized model for a given risk threshold <span class="math inline">\(\tau\)</span>. The final treatment recommendation is given by taking a majority vote from the individual models in the forest.</p>
<p>Below is the function call for constructing an rcRF model. In the next section we compare modeling results from rcDT and rcRF procedures graphically.</p>

<div class="sourceCode" id="cb10"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb10-1" title="1"><span class="kw">set.seed</span>(<span class="dv">2</span>)</a>
<a class="sourceLine" id="cb10-2" title="2">rcRF.fit &lt;-<span class="st"> </span><span class="kw">rcRF.select</span>(<span class="dt">data =</span> dat, </a>
<a class="sourceLine" id="cb10-3" title="3">                        <span class="dt">split.var =</span> <span class="dv">1</span><span class="op">:</span><span class="dv">10</span>, </a>
<a class="sourceLine" id="cb10-4" title="4">                        <span class="dt">efficacy =</span> <span class="st">&quot;y&quot;</span>, </a>
<a class="sourceLine" id="cb10-5" title="5">                        <span class="dt">risk =</span> <span class="st">&quot;r&quot;</span>,</a>
<a class="sourceLine" id="cb10-6" title="6">                        <span class="dt">col.trt =</span> <span class="st">&quot;trt&quot;</span>,</a>
<a class="sourceLine" id="cb10-7" title="7">                        <span class="dt">col.prtx =</span> <span class="st">&quot;prtx&quot;</span>,</a>
<a class="sourceLine" id="cb10-8" title="8">                        <span class="dt">risk.threshold =</span> <span class="fl">2.75</span>,</a>
<a class="sourceLine" id="cb10-9" title="9">                        <span class="dt">ntree =</span> <span class="dv">100</span>,</a>
<a class="sourceLine" id="cb10-10" title="10">                        <span class="dt">verbose =</span> <span class="ot">FALSE</span>)</a></code></pre></div>
</div>
<div id="predictions-for-rcdt-and-rcrf-models" class="section level2">
<h2>Predictions for rcDT and rcRF Models</h2>
<p>Prediction from rcDT and rcRF models can be obtained using the <code>predict.ITR()</code> function. The output from <code>predict.ITR()</code> is a list of several summaries from the model. If the model from which predictions are desired is an rcDT model, then of interest is <code>trt.pred</code> which gives the treatment recommendation vector for the prediction data. If the model is an rcRF, then <code>trt.pred</code> given the treatment recommendation based on majority voting. Also, <code>SummaryTreat</code> gives the probability of being recommended to active treatment (<span class="math inline">\(a = 1\)</span>) and <code>tree.votes</code> provides a matrix of votes from each tree.</p>

<div class="sourceCode" id="cb11"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb11-1" title="1">preds.rcDT &lt;-<span class="st"> </span><span class="kw">predict.ITR</span>(<span class="dt">fit =</span> rcDT.fit<span class="op">$</span>best.tree, </a>
<a class="sourceLine" id="cb11-2" title="2">                          <span class="dt">new.data =</span> dat, </a>
<a class="sourceLine" id="cb11-3" title="3">                          <span class="dt">split.var =</span> <span class="dv">1</span><span class="op">:</span><span class="dv">10</span>)<span class="op">$</span>trt.pred</a>
<a class="sourceLine" id="cb11-4" title="4">preds.rcRF &lt;-<span class="st"> </span><span class="kw">predict.ITR</span>(<span class="dt">fit =</span> rcRF.fit<span class="op">$</span>best.fit, </a>
<a class="sourceLine" id="cb11-5" title="5">                          <span class="dt">new.data =</span> dat, </a>
<a class="sourceLine" id="cb11-6" title="6">                          <span class="dt">split.var =</span> <span class="dv">1</span><span class="op">:</span><span class="dv">10</span>)<span class="op">$</span>trt.pred</a></code></pre></div>
<div class="figure">
<img src="data:image/png;base64,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" alt="Figure 2. Prediction Comparisons for rcDT and rcRF Models." />
<p class="caption">Figure 2. Prediction Comparisons for rcDT and rcRF Models.</p>
</div>
</div>
<div id="variable-importance" class="section level2">
<h2>Variable Importance</h2>
<p>Finally, we include a function to calculate the importance of a predictor relative to the total Lagrangian objective, the benefit scores, and the risk scores. This is accomplished by using the function <code>Variable.Importance.ITR()</code>. All three importance measures are calculated and returned. Here the <code>sort</code> argument is set to FALSE so that the importances are not individually ordered.</p>

<div class="sourceCode" id="cb12"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb12-1" title="1">VI &lt;-<span class="st"> </span><span class="kw">Variable.Importance.ITR</span>(<span class="dt">rcRF.fit =</span> rcRF.fit<span class="op">$</span>best.fit, </a>
<a class="sourceLine" id="cb12-2" title="2">                              <span class="dt">sort =</span> <span class="ot">FALSE</span>)</a>
<a class="sourceLine" id="cb12-3" title="3"><span class="kw">do.call</span>(cbind, VI)</a></code></pre></div>
<pre><code>##              VI VI.Efficacy     VI.Risk
## x1  0.285692680 0.347036411 0.030529479
## x2  0.663776169 0.598648287 0.899783370
## x3  0.007945324 0.008587843 0.011958512
## x4  0.007106902 0.007387884 0.012183190
## x5  0.005588163 0.005883192 0.007575996
## x6  0.010194106 0.011382871 0.010402283
## x7  0.006541436 0.007380822 0.009161867
## x8  0.006534948 0.006326655 0.010181776
## x9  0.004240365 0.004529210 0.005312606
## x10 0.002379908 0.002836826 0.002910921</code></pre>
<p>The variable importance results would be interpreted as <span class="math inline">\(X_1\)</span> and <span class="math inline">\(X_2\)</span> being important contributors to the overall rule (<code>VI</code>). Both <span class="math inline">\(X_1\)</span> and <span class="math inline">\(X_2\)</span> contribute similarly to efficacy (benefit) scores, but <span class="math inline">\(X_2\)</span> is strongly predictive of risk scores, i.e. important in keeping the risk constrained. </p>
</div>
<div id="references" class="section level2">
<h2>References</h2>
<p>[In Submission] Doubleday, K., Zhou, H., Fu, H., Zhou, J. (2020), “Risk Controlled Decision Trees and Random Forests for Precision Medicine”.</p>
</div>
</div>



<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>

</body>
</html>
