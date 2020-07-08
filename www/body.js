 $(document).ready(function(){

   $("#data_table .load-container").hide();
   $("#file_upload_btn").click(function(){
     $("#data_table .load-container").show();
   });
   
   $("#gobp_box .load-container").hide();
   $("#gsa_btn").click(function(){
     $("#gobp_box .load-container").show();
   });
   
    $("#gocc_box .load-container").hide();
   $("#gsa_btn").click(function(){
     $("#gocc_box .load-container").show();
   });
   
    $("#gomf_box .load-container").hide();
   $("#gsa_btn").click(function(){
     $("#gomf_box .load-container").show();
   });
   
    $("#kegg_box .load-container").hide();
   $("#gsa_btn").click(function(){
     $("#kegg_box .load-container").show();
   });
   

   content_width = $(".content").css('width');
   rightbar_width = $("#controlbar aside").css('width');
   size_table = parseInt(content_width) - parseInt(rightbar_width) - 30;
   size_plot = (size_table / 2) - 10;
      
   $(".nav a").click(function(){
      
      bodyClass = $('body').attr('class');
      if(bodyClass == 'skin-blue sidebar-mini control-sidebar-open'){
        
        dt_parent1 = $("#data_table").parent('div');
        dt_parent2 = dt_parent1.parent('div');
        dt_parent2.css('width',size_table);
        
        pt_parent1 = $("#plot_tabBox").parent('div');
        pt_parent2 = pt_parent1.parent('div');
        pt_parent2.css('width',size_table);
        
        //pca_plot_box
        pp_parent1 = $("#pca_plot_box").parent('div');
        pp_parent2 = pp_parent1.parent('div');
        pp_parent2.css('width',size_plot);
        
        //volcano_box
        vp_parent1 = $("#volcano_box").parent('div');
        vp_parent2 = vp_parent1.parent('div');
        vp_parent2.css('width',size_plot);
        
        //correlation_matrix_box
        cm_parent1 = $("#correlation_matrix_box").parent('div');
        cm_parent2 = cm_parent1.parent('div');
        cm_parent2.css('width',size_plot);
        
        //heatmap_box
        hm_parent1 = $("#heatmap_box").parent('div');
        hm_parent2 = hm_parent1.parent('div');
        hm_parent2.css('width',size_plot);
        
      } else{
         //data_table
        dt_parent1 = $("#data_table").parent('div');
        dt_parent2 = dt_parent1.parent('div');
        dt_parent2.css('width','');
        
          pt_parent1 = $("#plot_tabBox").parent('div');
        pt_parent2 = pt_parent1.parent('div');
        pt_parent2.css('width','');
        
        //pca_plot_box
        pp_parent1 = $("#pca_plot_box").parent('div');
        pp_parent2 = pp_parent1.parent('div');
        pp_parent2.css('width','');
        
        //volcano_box
        vp_parent1 = $("#volcano_box").parent('div');
        vp_parent2 = vp_parent1.parent('div');
        vp_parent2.css('width','');
        
        //correlation_matrix_box
        cm_parent1 = $("#correlation_matrix_box").parent('div');
        cm_parent2 = cm_parent1.parent('div');
        cm_parent2.css('width','');
        
        //heatmap_box
        hm_parent1 = $("#heatmap_box").parent('div');
        hm_parent2 = hm_parent1.parent('div');
        hm_parent2.css('width','');
      }
      
   });
 });