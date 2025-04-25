#' SimpleHeatmap
#'
#' @param Cindex_mat A matrix containing C-index values for visualization. Rows correspond to items, and columns correspond to cohorts.
#' @param avg_Cindex A numeric vector representing the average C-index values for rows in Cindex_mat.
#' @param CohortCol A named vector or list specifying the colors for each cohort in Cindex_mat. Names must match the column names of Cindex_mat.
#' @param barCol A vector specifying the fill colors for the bar plot representing avg_Cindex. Each element corresponds to a row in avg_Cindex.
#' @param cellwidth A numeric value specifying the width of each cell in the heatmap (default: 1).
#' @param cellheight A numeric value specifying the height of each cell in the heatmap (default: 0.5).
#' @param cluster_columns A logical value or clustering method specifying whether to cluster the columns of the heatmap. Defaults to user-specified input.
#' @param cluster_rows A logical value or clustering method specifying whether to cluster the rows of the heatmap. Defaults to user-specified input.
#'
#' @return A heatmap object visualizing the C-index matrix, annotated with cohort labels and bar plots for row-level averages.
#' @export
#'
#' @examples
SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)

  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)

  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha,
          top_annotation = col_ha,
          # col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 黄绿配色
          col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
          rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, # 不进行聚类，无意义
          show_column_names = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}



