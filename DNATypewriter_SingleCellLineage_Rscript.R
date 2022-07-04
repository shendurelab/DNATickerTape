## 2022.07.04 Plotting cell lineage tree from the input of cell-by-edit csv files.


# Load R packages
library(ape)
library(phangorn)
library(tidyverse)
library(ggdendro)
library(dendextend)

# =====================================================================
# Generating the wider cell-by-59EditSites file from cell-by-5EditSites
# =====================================================================


edit_table_by_5 = read.csv("Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))


### Filtering TargetBC set based on the number of retained cell (13 TargetBC, 3257 cells)
TargetBC_freq <- data.frame(table(edit_table_by_5$TargetBC)) %>% arrange(desc(Freq))
retained_cells <- rep(0,17)
for (ii in 1:17){
  current_TargetBC_set <- as.character(TargetBC_freq$Var1)[1:ii]
  Cell_filter_by_TargetBC <- select(edit_table_by_5, c('Cell','TargetBC')) %>%
    filter(TargetBC %in% current_TargetBC_set) %>%
    mutate(count = 1) %>%
    pivot_wider(names_from = TargetBC, values_from = 'count', values_fill = 0)
  retained_cells[ii] <- sum(rowSums(Cell_filter_by_TargetBC[,2:(ii+1)]) == ii)
}
retained_cells <- data.frame(retained_cells)
retained_cells$numTargetBC <- 1:17
retained_cells$retained_cells <- as.numeric(as.character(retained_cells$retained_cells))
targets_to_use <- as.character(TargetBC_freq$Var1[1:13])


### Generating the cell_list for 3,257 cells
Cell_filter_13set <- select(edit_table_by_5, c('Cell','TargetBC')) %>%
  filter(TargetBC %in% targets_to_use) %>%
  mutate(count = 1) %>%
  pivot_wider(id_cols = 'Cell',names_from = TargetBC, values_from = 'count', values_fill = 0)
Cell_filter_13set$sums <- rowSums(Cell_filter_13set[,2:14])
Cell_filter_13set <- filter(Cell_filter_13set, sums == 13)
cell_list_3257 <- Cell_filter_13set$Cell
edit_table_3257 <- filter(edit_table_by_5, TargetBC %in% targets_to_use & Cell %in% cell_list_3257)

# Marking unedited sites with 'None', and non-existing sites (contracted to 4xTAPE or 2xTAPE) with NA
edit_table_3257[is.na(edit_table_3257)] <- 'None'
edit_table_3257[edit_table_3257$TargetBC == 'TGGACGAC',7] <- NA
edit_table_3257[edit_table_3257$TargetBC == 'TTTCGTGA',7] <- NA
edit_table_3257[edit_table_3257$TargetBC == 'TGGTTTTG',7] <- NA
edit_table_3257[edit_table_3257$TargetBC == 'TTCACGTA',5:7] <- NA


# edit_cell_table_65 = Ordered cell-by-65EditSites table, including non-existing sites before contraction
edit_cell_table_65 <- select(edit_table_3257, -nUMI) %>%
  pivot_longer(cols = c('Site1','Site2','Site3','Site4','Site5'), names_to = 'Sites', values_to ='Insert') %>%
  pivot_wider(id_cols = Cell, names_from = c(TargetBC,Sites), names_sep = ".", values_from = Insert)
edit_cell_table_65 <- arrange(edit_cell_table_65,Cell) %>%
  select(order(colnames(edit_cell_table_65)))


sub_edit65 <- as.matrix(select(edit_cell_table_65,-Cell))
sub_edit65[is.na(sub_edit65)] <- 'None'
rownames(sub_edit65) <- edit_cell_table_65$Cell

sub_edit59 <- as.matrix(select(edit_cell_table_65,-c('Cell','TGGACGAC.Site5','TTTCGTGA.Site5','TGGTTTTG.Site5',
                                              'TTCACGTA.Site3','TTCACGTA.Site4','TTCACGTA.Site5')))
rownames(sub_edit59) <- edit_cell_table_65$Cell
cell_list <- edit_cell_table_65$Cell


# =====================================================================
# Generating the phylogenetic tree based on edits
# =====================================================================

# shared_edit_matrix = Counting all shared edits per cell-pair, consistent with the sequential editing on DNA Tape



# Function for calculating shared_edit_matrix
fun_shared_edit_matrix <- function(x) {
  #if (1 == 1){
  sub_edit65 <- x
  sub_edit65[sub_edit65 == 'None'] <- 1:filter(as.data.frame(table(sub_edit65)), sub_edit65 == 'None')$Freq
  ncell <- nrow(sub_edit65)
  cell_list <- sort(rownames(sub_edit65))
  shared_edit_matrix <- matrix(0, ncell,ncell)
  colnames(shared_edit_matrix) <- cell_list
  rownames(shared_edit_matrix) <- cell_list
  for (ii in 1:(ncell)){
    cell1 <- cell_list[ii]
    for (jj in (ii):ncell){
      cell2 <- cell_list[jj]
      for (kk in seq(0,(dim(sub_edit65)[2]-5),5)){
        if (sub_edit65[cell1,(kk+1)] == sub_edit65[cell2,(kk+1)]){
          shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
          if (sub_edit65[cell1,(kk+2)] == sub_edit65[cell2,(kk+2)]){
            shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
            if (sub_edit65[cell1,(kk+3)] == sub_edit65[cell2,(kk+3)]){
              shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
              if (sub_edit65[cell1,(kk+4)] == sub_edit65[cell2,(kk+4)]){
                shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                if (sub_edit65[cell1,(kk+5)] == sub_edit65[cell2,(kk+5)]){
                  shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                }
              }
            }
          }
        }
      }
      shared_edit_matrix[jj,ii] <- shared_edit_matrix[ii,jj]
    }
  }
  return(shared_edit_matrix)
}

#shared_edit_matrix <- fun_shared_edit_matrix(sub_edit65) #15-30 min on computer; once done, saved and loaded for the future use
shared_edit_matrix <- read.csv('shared_edit_matrix_3257.csv', stringsAsFactors = F, header = T)
shared_edit_matrix <- as.matrix(shared_edit_matrix)
diag(shared_edit_matrix) <- 59

distance_matrix <- 59 - shared_edit_matrix # Phylogenetic distance caludated as (# of possible sites - # of shared sites)
distance_matrix <- as.matrix(distance_matrix)
tree <- as.phylo(hclust(as.dist(distance_matrix), "average")) # tree built using UPGMA


# Bootstrapping = takes ~1 day. Once done, saved and loaded for the future use
#bstrees_3257 <- boot.phylo(tree, sub_edit65, fun_order_phylo,  mc.cores = 7, trees = TRUE, block = 5)$trees
#bootstrap_value_3257 <- prop.clades(tree, bstrees_3257, rooted = TRUE)
#write.table(data.frame(bootstrap_value_3257),'TreeBSvalues_calc_3257.csv',sep=',')
bootstrap_value_3257 <- read.csv('TreeBSvalues_calc_3257.csv', stringsAsFactors = F, header = T)
bootstrap_value_3257 <- bootstrap_value_3257$clad_3257



# =====================================================================
# Plotting the entire tree
# =====================================================================


### Plot size
# load color scheme for each edits
color_legend <- read.csv('edit_color_legend.csv', stringsAsFactors = F, header = T)
cols<-color_legend$cols
names(cols) <- color_legend$label
cell_list <- edit_cell_table_65$Cell
plot_height <- 3257
plot_width <- 50




# Calculate the node position for writing clad info
nd <- node.depth(tree, method = 1)
nde <- node.depth.edgelength(tree)
nh <- node.height(tree)
nodes_positions <- bind_cols(as.data.frame(nd),as.data.frame(nde),as.data.frame(nh))
nodes_positions$x <- nodes_positions$nde / max(nodes_positions$nde) * plot_width
nodes_positions$y <- (nodes_positions$nh- min(nodes_positions$nh)) / (max(nodes_positions$nh) - min(nodes_positions$nh)) * plot_height
nodes_positions <- filter(nodes_positions, nd != 1)
nodes_positions <- bind_cols(nodes_positions, as.data.frame(bootstrap_value_3257))


### Get coordinates for the dendogram segments
tree_dend <- as.dendrogram(hclust(as.dist(distance_matrix), "average"))
dend_data <- dendro_data(tree_dend, type = "rectangle")

# Rescale X and Y
branch_positions <- dend_data$segments

x_min <- min(branch_positions$x, branch_positions$xend)
x_max <- max(branch_positions$x, branch_positions$xend)
y_min <- min(branch_positions$y, branch_positions$yend)
y_max <- max(branch_positions$y, branch_positions$yend)
branch_positions$x <- (branch_positions$x - x_min)/(x_max - x_min) * plot_height
branch_positions$xend <- (branch_positions$xend - x_min)/(x_max - x_min) * plot_height
branch_positions$y <- (1-(branch_positions$y - y_min)/(y_max - y_min) )* plot_width
branch_positions$yend <- (1-(branch_positions$yend - y_min)/(y_max - y_min) )* plot_width
colnames(branch_positions) <- c('y','x','yend','xend')

# Assign label positions
branch_labels <- dend_data$labels
label_pos <- filter(branch_positions, xend == plot_width) %>%
  arrange(yend)
branch_labels$x <- label_pos$xend
branch_labels$y <- label_pos$yend

cell_list <- branch_labels$label
label_pos <-select(branch_labels, c('label','y'))
colnames(label_pos) <- c('cell_name','ylabel')



#Create branch_label_ends file
branch_label_ends <- select(branch_labels, -x)
colnames(branch_label_ends) = c('y','cell_name')
branch_label_ends$x1 <- 0
branch_label_ends$x2 <- 1


# Making cigar file
edit_cell_table_65 <- arrange(edit_cell_table_65,match(Cell, branch_label_ends$cell_name))
sub_edit_table <- as.matrix(select(edit_cell_table_65,-Cell))
rownames(sub_edit_table) <- edit_cell_table_65$Cell
sub_edit_table[is.na(sub_edit_table)] <- 'None'


branch_label_ends$InsertBC <- sub_edit_table[,1]

# Make the edit_plot_whole file 
shift_x_width = 1
edit_plot_whole = branch_label_ends
terminal_branch_info = edit_plot_whole
short_TAPE_list <- c(40,45,53,54,55,65) #contracted DNA tape locations

for (i in 2:ncol(sub_edit_table)){
  temp_dataframe = terminal_branch_info
  if (i %% 5 == 1){
    temp_dataframe$x1 = temp_dataframe$x1 + shift_x_width
    temp_dataframe$x2 = temp_dataframe$x2 + shift_x_width
  }
  if (i %in% short_TAPE_list){
    terminal_branch_info = temp_dataframe
  }else{
    temp_dataframe$x1 = temp_dataframe$x1 + shift_x_width
    temp_dataframe$x2 = temp_dataframe$x2 + shift_x_width
    temp_dataframe$InsertBC = sub_edit_table[,i]
    terminal_branch_info = temp_dataframe
    edit_plot_whole = rbind(edit_plot_whole,temp_dataframe)
    #print(nrow(edit_plot_whole))
  }
  
}
edit_plot_whole_sub <- edit_plot_whole
edit_plot_whole_sub <- left_join(edit_plot_whole_sub, label_pos, by = 'cell_name')
edit_plot_whole_sub$x1 <- edit_plot_whole_sub$x1 + max(branch_positions$xend)
edit_plot_whole_sub$x2 <- edit_plot_whole_sub$x2 + max(branch_positions$xend)
branch_labels$xend <- max(edit_plot_whole_sub$x2)


##Add target numbers:
square_width = edit_plot_whole_sub[1,]$x2 - edit_plot_whole_sub[1,]$x1
TargetBC_label = colnames(sub_edit59)
xpos_of_TargetBC_label = unique(edit_plot_whole_sub$x2)-(.5*square_width)
TargetBC_position_table = as.data.frame(cbind(TargetBC_label,xpos_of_TargetBC_label))
TargetBC_position_table$xpos_of_TargetBC_label = as.numeric(TargetBC_position_table$xpos_of_TargetBC_label)
TargetBC_position_table$y_pos = min(edit_plot_whole_sub$ylabel) - 1



tree_plot = ggplot() + geom_point()
tree_plot = tree_plot + geom_segment(data = branch_positions, aes(x = x, y = y, xend = xend, yend = yend))
tree_plot = tree_plot + 
  theme_minimal() + theme(legend.position="none") + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank()) +
  theme(panel.background = element_rect(fill = 'white', color = 'white')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tree_plot = tree_plot + geom_rect(data =edit_plot_whole_sub, mapping=aes(xmin = x2, xmax = x2+1, ymin=ylabel-.4, ymax=ylabel+.4, fill = factor(InsertBC)), color = "black", size = 0.2)
tree_plot = tree_plot + scale_fill_manual(values = c(cols))
tree_plot = tree_plot + geom_text(aes(x = branch_labels$x, y = branch_labels$y, label = round(branch_labels$y, digits = 0)), hjust = 0, nudge_x = 0, size = 1)
tree_plot = tree_plot + geom_text(aes(x = branch_labels$xend+2, y = branch_labels$y, label = branch_labels$label), hjust = 0, nudge_x = 0, size = 1)
tree_plot = tree_plot + geom_text(aes(x = nodes_positions$x, y = nodes_positions$y, label = nodes_positions$bootstrap_value_3257), color = "blue", hjust = 0, nudge_x = 0, size = 1.2)
tree_plot = tree_plot + geom_text(aes(x = TargetBC_position_table$xpos_of_TargetBC_label+1, y = TargetBC_position_table$y_pos-7, label = TargetBC_position_table$TargetBC_label), hjust = 0, nudge_x = 0, nudge_y = -.3, size = 1.2, angle = 90)
tree_plot = tree_plot + xlim(0, max(edit_plot_whole_sub$x2 + 15))

ggsave(filename = 'Supplementary_File_1_tree3257.pdf', plot = tree_plot, width =6, height = 160, limitsize = FALSE)
