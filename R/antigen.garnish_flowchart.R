

antigen.garnish_flowchart <<- try({

glue::glue("digraph boxes_and_circles {{

graph [overlap = true, bgcolor = transparent, rankdir = LR]

node [margin = 0.2, shape = box, fontsize = 18, fontname = '{Font}', fontcolor = '#000000', color = '#000000', fillcolor = '#82b1ff:#8c9eff', style = radial]

  sv [label = 'Somatic variants\n(VCF file)']
  se [label = 'Annotate (SnpEff)']
  epi [label = 'Predicted neoepitopes']

node [margin = 0.2, shape = box, fontname = '{Font}', fontcolor = '#000000', fontsize = 18, color = '#000000', fillcolor = '#ffd180:#ffab40', style = radial]

  ag [label = 'antigen.garnish']

sv -> se -> ag -> epi

edge [penwidth = 1.5, color = '#000000', arrowsize = 1]

}}")}, silent = TRUE)
