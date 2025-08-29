library(DiagrammeR)
library(glue)
library(readxl)
library(DiagrammeRsvg)
library(rsvg)
library(tidyverse)

arq <- 'tr_t_veh'

file_path <- glue::glue("C:/Users/AsRock/Documents/{arq}.xlsx")

tr <- data.frame(readxl::read_excel(file_path))

tr[, 1:(ncol(tr) - 1)] <- lapply(tr[, 1:(ncol(tr) - 1)], function(col) {
  ifelse(col < 5, 0, col)
})

tr[, 1:(ncol(tr) - 1)] <- lapply(tr[, 1:(ncol(tr) - 1)], function(x) x*0.5)

tr[, ncol(tr)] <- log(tr[, ncol(tr)])

tr <- lapply(tr, function(x) {
  ifelse(x < 0, 0, x)
})

tr <- data.frame(tr)

rownames(tr) <- colnames(tr[1:ncol(tr)-1])

c <- tr['c', ncol(tr)]
g <- tr['g', ncol(tr)]
i <- tr['i', ncol(tr)]
l <- tr['l', ncol(tr)]
r <- tr['r', ncol(tr)]
t <- tr['t', ncol(tr)]
w <- tr['w', ncol(tr)]
b <- tr['b', ncol(tr)]
o <- tr['o', ncol(tr)]

wc = '#29b6f6'
tc = '#3ad740'
cc = '#b6b70f'
gc = '#9c27b0'
ic = '#3949ab'
rc = '#f48fb1'
lc = '#f44336'
oc = '#cacfd2'
bc = '#f4d03f' 



eto <- DiagrammeR::grViz(glue::glue(
  .open = "{{",
  .close = "}}", 
  "digraph {

  graph [rankdir = TB,
         layout = circo, 
         fontname = Arial,
         splines = false,
         overlap = false,
         ranksep = 1,
         height = 5,
         width = 0.5,
         normalize = 0
         ]
  
  node [shape = circle, 
        fontsize = 50,
        style = filled,
        label = ''
        ]
       
	edge [arrowhead = ldiamond, 
	      arrowsize = 0.1,
	      headclip = true,
	      len = 6,
	      ]
	      
   
w [xlabel = 'W', style='filled' fillcolor= '{{wc}}', width = {{w}}]
t [xlabel = 'T', style='filled' fillcolor= '{{tc}}', width = {{t}}]
c [xlabel = 'C', style='filled' fillcolor= '{{cc}}', width = {{c}}]
g [xlabel = 'G', style='filled' fillcolor= '{{gc}}', width = {{g}}]
i [xlabel = 'I', style='filled' fillcolor= '{{ic}}', width = {{i}}]
r [xlabel = 'R', style='filled' fillcolor= '{{rc}}', width = {{r}}]
l [xlabel = 'L', style='filled' fillcolor= '{{lc}}', width = {{l}}]
o [xlabel = 'O', style='filled' fillcolor= '{{oc}}', width = {{o}}]
b [xlabel = 'B', style='filled' fillcolor= '{{bc}}', width = {{b}}]  

g -> g [color= '{{gc}}', penwidth = {{tr['g','g']}}, weight = {{g}}]
g -> l [color= '{{gc}}', penwidth = {{tr['g','l']}}, weight = {{g}}]
g -> r [color= '{{gc}}', penwidth = {{tr['g','r']}}, weight = {{g}}]
g -> t [color= '{{gc}}', penwidth = {{tr['g','t']}}, weight = {{g}}]
g -> w [color= '{{gc}}', penwidth = {{tr['g','w']}}, weight = {{g}}]
g -> i [color= '{{gc}}', penwidth = {{tr['g','i']}}, weight = {{g}}]
g -> o [color= '{{gc}}', penwidth = {{tr['g','o']}}, weight = {{g}}]
g -> c [color= '{{gc}}', penwidth = {{tr['g','c']}}, weight = {{g}}]
g -> b [color= '{{gc}}', penwidth = {{tr['g','b']}}, weight = {{g}}]

l -> g [color= '{{lc}}', penwidth = {{tr['l','g']}}, weight = {{l}}]
l -> l [color= '{{lc}}', penwidth = {{tr['l','l']}}, weight = {{l}}]
l -> r [color= '{{lc}}', penwidth = {{tr['l','r']}}, weight = {{l}}]
l -> t [color= '{{lc}}', penwidth = {{tr['l','t']}}, weight = {{l}}]
l -> w [color= '{{lc}}', penwidth = {{tr['l','w']}}, weight = {{l}}]
l -> i [color= '{{lc}}', penwidth = {{tr['l','i']}}, weight = {{l}}]
l -> o [color= '{{lc}}', penwidth = {{tr['l','o']}}, weight = {{l}}]
l -> c [color= '{{lc}}', penwidth = {{tr['l','c']}}, weight = {{l}}]
l -> b [color= '{{lc}}', penwidth = {{tr['l','b']}}, weight = {{l}}]

r -> g [color= '{{rc}}', penwidth = {{tr['r','g']}}, weight = {{r}}]
r -> l [color= '{{rc}}', penwidth = {{tr['r','l']}}, weight = {{r}}]
r -> r [color= '{{rc}}', penwidth = {{tr['r','r']}}, weight = {{r}}]
r -> t [color= '{{rc}}', penwidth = {{tr['r','t']}}, weight = {{r}}]
r -> w [color= '{{rc}}', penwidth = {{tr['r','w']}}, weight = {{r}}]
r -> i [color= '{{rc}}', penwidth = {{tr['r','i']}}, weight = {{r}}]
r -> o [color= '{{rc}}', penwidth = {{tr['r','o']}}, weight = {{r}}]
r -> c [color= '{{rc}}', penwidth = {{tr['r','c']}}, weight = {{r}}]
r -> b [color= '{{rc}}', penwidth = {{tr['r','b']}}, weight = {{r}}]

t -> g [color= '{{tc}}', penwidth = {{tr['t','g']}}, weight = {{t}}]
t -> l [color= '{{tc}}', penwidth = {{tr['t','l']}}, weight = {{t}}]
t -> r [color= '{{tc}}', penwidth = {{tr['t','r']}}, weight = {{t}}]
t -> t [color= '{{tc}}', penwidth = {{tr['t','t']}}, weight = {{t}}]
t -> w [color= '{{tc}}', penwidth = {{tr['t','w']}}, weight = {{t}}]
t -> i [color= '{{tc}}', penwidth = {{tr['t','i']}}, weight = {{t}}]
t -> o [color= '{{tc}}', penwidth = {{tr['t','o']}}, weight = {{t}}]
t -> c [color= '{{tc}}', penwidth = {{tr['t','c']}}, weight = {{t}}]
t -> b [color= '{{tc}}', penwidth = {{tr['t','b']}}, weight = {{t}}]

w -> g [color= '{{wc}}', penwidth = {{tr['w','g']}}, weight = {{w}}]
w -> l [color= '{{wc}}', penwidth = {{tr['w','l']}}, weight = {{w}}]
w -> r [color= '{{wc}}', penwidth = {{tr['w','r']}}, weight = {{w}}]
w -> t [color= '{{wc}}', penwidth = {{tr['w','t']}}, weight = {{w}}]
w -> w [color= '{{wc}}', penwidth = {{tr['w','w']}}, weight = {{w}}]
w -> i [color= '{{wc}}', penwidth = {{tr['w','i']}}, weight = {{w}}]
w -> o [color= '{{wc}}', penwidth = {{tr['w','o']}}, weight = {{w}}]
w -> c [color= '{{wc}}', penwidth = {{tr['w','c']}}, weight = {{w}}]
w -> b [color= '{{wc}}', penwidth = {{tr['w','b']}}, weight = {{w}}]

i -> g [color= '{{ic}}', penwidth = {{tr['i','g']}}, weight = {{i}}]
i -> l [color= '{{ic}}', penwidth = {{tr['i','l']}}, weight = {{i}}]
i -> r [color= '{{ic}}', penwidth = {{tr['i','r']}}, weight = {{i}}]
i -> t [color= '{{ic}}', penwidth = {{tr['i','t']}}, weight = {{i}}]
i -> w [color= '{{ic}}', penwidth = {{tr['i','w']}}, weight = {{i}}]
i -> i [color= '{{ic}}', penwidth = {{tr['i','i']}}, weight = {{i}}]
i -> o [color= '{{ic}}', penwidth = {{tr['i','o']}}, weight = {{i}}]
i -> c [color= '{{ic}}', penwidth = {{tr['i','c']}}, weight = {{i}}]
i -> b [color= '{{ic}}', penwidth = {{tr['i','b']}}, weight = {{i}}]

c -> g [color= '{{cc}}', penwidth = {{tr['c','g']}}, weight = {{c}}]
c -> l [color= '{{cc}}', penwidth = {{tr['c','l']}}, weight = {{c}}]
c -> r [color= '{{cc}}', penwidth = {{tr['c','r']}}, weight = {{c}}]
c -> t [color= '{{cc}}', penwidth = {{tr['c','t']}}, weight = {{c}}]
c -> w [color= '{{cc}}', penwidth = {{tr['c','w']}}, weight = {{c}}]
c -> i [color= '{{cc}}', penwidth = {{tr['c','i']}}, weight = {{c}}]
c -> o [color= '{{cc}}', penwidth = {{tr['c','o']}}, weight = {{c}}]
c -> c [color= '{{cc}}', penwidth = {{tr['c','c']}}, weight = {{c}}]
c -> b [color= '{{cc}}', penwidth = {{tr['c','b']}}, weight = {{c}}]

o -> g [color= '{{oc}}', penwidth = {{tr['o','g']}}, weight = {{o}}]
o -> l [color= '{{oc}}', penwidth = {{tr['o','l']}}, weight = {{o}}]
o -> r [color= '{{oc}}', penwidth = {{tr['o','r']}}, weight = {{o}}]
o -> t [color= '{{oc}}', penwidth = {{tr['o','t']}}, weight = {{o}}]
o -> w [color= '{{oc}}', penwidth = {{tr['o','w']}}, weight = {{o}}]
o -> i [color= '{{oc}}', penwidth = {{tr['o','i']}}, weight = {{o}}]
o -> o [color= '{{oc}}', penwidth = {{tr['o','o']}}, weight = {{o}}]
o -> c [color= '{{oc}}', penwidth = {{tr['o','c']}}, weight = {{o}}]
o -> b [color= '{{oc}}', penwidth = {{tr['o','b']}}, weight = {{o}}]

b -> g [color= '{{bc}}', penwidth = {{tr['b','g']}}, weight = {{b}}]
b -> l [color= '{{bc}}', penwidth = {{tr['b','l']}}, weight = {{b}}]
b -> r [color= '{{bc}}', penwidth = {{tr['b','r']}}, weight = {{b}}]
b -> t [color= '{{bc}}', penwidth = {{tr['b','t']}}, weight = {{b}}]
b -> w [color= '{{bc}}', penwidth = {{tr['b','w']}}, weight = {{b}}]
b -> i [color= '{{bc}}', penwidth = {{tr['b','i']}}, weight = {{b}}]
b -> o [color= '{{bc}}', penwidth = {{tr['b','o']}}, weight = {{b}}]
b -> c [color= '{{bc}}', penwidth = {{tr['b','c']}}, weight = {{b}}]
b -> b [color= '{{bc}}', penwidth = {{tr['b','b']}}, weight = {{b}}]


  }"))

eto

eto %>%
  export_svg() %>%
  charToRaw %>% 
  rsvg_pdf(glue("{arq}.svg"))
