$(document).ready(function(){
  
  //$(window).resize(function(e) {
    //dimension[0] = window.innerWidth;
    //dimension[1] = window.innerHeight;
   //Shiny.onInputChange("dimension", dimension);
  //});
  
  var data_table = "hide";
  $("#data_table .load-container").hide();
  
  $("#file_upload_btn").click(function(){
    if(data_table == "hide"){
     $("#data_table .load-container").show();
     data_table = "show";
    } else{
      $("#uploaded_file_header").attr('style','visibility:hidden');
    }
  });
  
  $(".btn-file").click(function(){
    if(data_table != "hide"){
      $("#uploaded_file_header").attr('style','visibility:hidden');
    }
  });

  $("#gobp_gsa_box .load-container").hide();
  $("#gocc_gsa_box .load-container").hide();
  $("#gomf_gsa_box .load-container").hide();
  $("#kegg_gsa_box .load-container").hide();
  
  $("#gsa_btn").click(function(){
    $("#gobp_gsa_box .load-container").show();
  });
  
  $("#gsa_btn").click(function(){
    $("#gocc_gsa_box .load-container").show();
  });
  
  $("#gsa_btn").click(function(){
    $("#gomf_gsa_box .load-container").show();
  });
  
  $("#gsa_btn").click(function(){
    $("#kegg_gsa_box .load-container").show(); 
  });
  
  $("#gobp_gsea_box .load-container").hide();
  $("#gocc_gsea_box .load-container").hide();
  $("#gomf_gsea_box .load-container").hide();
  $("#kegg_gsea_box .load-container").hide();
  
  $("#gsea_btn").click(function(){
    $("#gobp_gsea_box .load-container").show();
  });
   
  $("#gsea_btn").click(function(){
    $("#gocc_gsea_box .load-container").show();
  });
   
  $("#gsea_btn").click(function(){
    $("#gomf_gsea_box .load-container").show();
  });
   
   $("#gsea_btn").click(function(){
     $("#kegg_gsea_box .load-container").show();
   });
  
  $("#ppi_box .load-container").hide();
  $("#ppi_btn").click(function(){
     $("#ppi_box .load-container").show();
   });
  
  var content_width = $(".content").css('width');
  var rightbar_width = $("#controlbar aside").css('width');
  var size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
  var size_plot = (size_table / 2) - 10;
 
 
  $( window ).resize( function() {
    content_width = $(".content").css('width');
    rightbar_width = $("#controlbar aside").css('width');
    size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
    size_plot = (size_table / 2) - 10;

  });

  var dimension = [0, 0];
  $("#zoom_pathway_btn").click(function(){
    dimension[0] = 1344;
    dimension[1] = 756;
      Shiny.onInputChange("dimension", dimension);
  });
 
  $(".nav.navbar-nav a").click(function(){
    content_width = $(".content").css('width');
    rightbar_width = $("#controlbar aside").css('width');
    size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
    size_plot = (size_table / 2) - 10;
    
    bodyClass = $('body').attr('class');
    if(bodyClass == 'skin-blue sidebar-mini control-sidebar-open'){
      
      dt_parent1 = $("#data_table").parent('div');
      dt_parent2 = dt_parent1.parent('div');
      dt_parent2.css('width',size_table);
      
      pt_parent1 = $("#plot_tabBox").parent('div');
      pt_parent1.css('height',$("#plot_tabBox").css('height'));
      pt_parent2 = pt_parent1.parent('div');
      pt_parent2.css('width',size_table);
      
      //pca_plot_box
      pp_parent1 = $("#pca_plot_box").parent('div');
      pp_parent1.css('height',$("#pca_plot_box").css('height'));
      pp_parent2 = pp_parent1.parent('div');
      pp_parent2.css('width',size_plot);
      
      //volcano_box
      vp_parent1 = $("#volcano_box").parent('div');
      var vp_height = $("#volcano_box").css('height')+$("#volcano_info").css('height')+$(".download_volcano").css('height');
      vp_parent1.css('height',vp_height);
      vp_parent2 = vp_parent1.parent('div');
      vp_parent2.css('width',size_plot);
      
      //correlation_matrix_box
      cm_parent1 = $("#correlation_matrix_box").parent('div');
      cm_parent1.css('height',$("#correlation_matrix_box").css('height'));
      cm_parent2 = cm_parent1.parent('div');
      cm_parent2.css('width',size_plot);
      
      //heatmap_box
      hm_parent1 = $("#heatmap_box").parent('div');
      hm_parent1.css('height',$("#heatmap_box").css('height'));
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