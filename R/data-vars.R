#' @title Coding Class
#' 
#' @description Coding classes as a vector of characters.
#' 
#' @details This object is called by [GAMBLR.data::get_coding_ssm] and [GAMBLR.data::assign_cn_to_ssm].
#' 
#' @noRd
#' 
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                 "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", 
                 "Nonstop_Mutation", "Silent", "Splice_Region", 
                 "Splice_Site", "Targeted_Region", "Translation_Start_Site")
