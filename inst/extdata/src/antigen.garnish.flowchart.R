   antigen.garnish_flowchart <- glue::glue("digraph boxes_and_circles {{

          graph [overlap = true, bgcolor = transparent, rankdir = LR]

          node [margin = 0.2, shape = box, fontsize = 18, fontname = '{Font}', fontcolor = '#000000', color = '#000000', fillcolor = '#82b1ff:#8c9eff', style = radial]

            input [label = 'Somatic variants\n(VCF file)\nor peptides']
            epi [label = 'Top predicted\nneoepitopes']

          node [margin = 0.2, face = bold, shape = box, fontname = '{Font}', fontcolor = '#000000', fontsize = 22, color = '#000000', fillcolor = '#ffd180:#ff9100', style = radial]

            ag [label = 'antigen.garnish']

          input -> ag -> epi

          edge [penwidth = 1.5, color = '#000000', arrowsize = 1]

          }}")

        antigen.garnish_flowchart %<>% DiagrammeR::grViz(.)
        write(DiagrammeRsvg::export_svg(antigen.garnish_flowchart),
          file = paste0(
            DROPBOX,
            "S3/antigen.garnish_flowchart.svg"
          ),
          ncolumns = 1,
          append = FALSE,
          sep = " "
        )
