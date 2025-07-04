<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="icon" href="https://raw.githubusercontent.com/artic-network/website/master/assets/images/artic-logo-small.svg">

    <title>squirrel</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha384-HSMxcRTRxnN+Bdg0JdbxYKrThecOKuH5zCYotlSAcp1+c8xmyTe9GYg1l9a69psu" crossorigin="anonymous">
    <script src="https://code.jquery.com/jquery-1.12.4.min.js" integrity="sha384-nvAa0+6Qg9clwYCGGPpDQLVpLNn0fRaROjHqs13t4Ggj3Ez50XnGQqc/r8MhnRDZ" crossorigin="anonymous"></script>
    <link href="https://cdn.datatables.net/1.10.25/css/jquery.dataTables.min.css" rel="stylesheet" type="text/css" />
    <script src="https://cdn.datatables.net/1.10.25/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/dataTables.buttons.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.html5.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.print.min.js"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha384-aJ21OjlMXNL5UyIl/XNwTMqvzeRMZH2w8c5cRVpzpU8Y5bApTppSuUkhZXN0VxHd" crossorigin="anonymous"></script>
    <script src="https://cdn.jsdelivr.net/gh/rambaut/figtree.js@9880/dist/figtree.umd.js"></script>
    <script src="https://d3js.org/d3.v6.min.js"></script>
    <script src="https://sharonchoong.github.io/svg-exportJS/svg-export.min.js"></script>
    <script src="https://unpkg.com/canvg/lib/umd.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/pdfkit/js/pdfkit.min.js"></script>
    <script src="https://github.com/devongovett/blob-stream/releases/download/v0.1.3/blob-stream.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/svg-to-pdfkit/source.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega@5.16.0"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@4.15.0"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@6.11.1"></script>

    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <% themeColor = "#419595" %>
    <style>
      body {
        padding-top: 50px;
        font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif;
      }
      table text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
      }
      header {
          display: block;
          text-align: right;
      
      }
      svg {width: 90%; height: auto;}
      .row {
          display: flex;
        }
        .column {
          padding: 10px;
          flex: 50%;
        }
      .accordion {
          background-color: #eee;
          color: #444;
          cursor: pointer;
          padding: 13px;
          width: 100%;
          border: none;
          text-align: left;
          outline: none;
          transition: 0.4s;
        }

        .active, .accordion:hover {
          background-color: ${themeColor};
          color: white;
        }

        .accordion:after {
          content: '\002B';
          color: white;
          font-weight: bold;
          float: right;
          margin-left: 5px;
        }

        .active:after {
          content: "\2212";
        }

        .panel {
          padding: 0 13px;
          background-color: white;
          max-height: 0;
          overflow: hidden;
          transition: max-height 0.2s ease-out;
        }
      .center {
          display: block;
          margin-left: auto;
          margin-right: auto;
          width: 50%;
          }
      .node-background{
          fill:dimgrey;
          stroke:dimgrey;
      }
      .node circle{
        stroke-width:0;
        cursor:pointer;
        /* fill:#7178bc; */
        stroke:dimgrey;
        }
      .node circle.selected{
        stroke-width:0;
        cursor:pointer;
        fill:${themeColor};
        stroke:dimgrey;
        }
      .node-background.query_boolean-True{
          stroke:${themeColor};
      }
      .node.query_boolean-True circle{
        stroke:${themeColor};
      }
      .node.query_boolean-True circle.selected{
        stroke:${themeColor};
      }
      .node-background.query_boolean-True circle.selected{
          stroke:${themeColor};
      }
      .node.query_boolean-True.hovered circle{
          stroke:${themeColor};
      }
      .node rect{
        stroke-width:2;
        fill:${themeColor};
        stroke:dimgrey;
      }
      .svg-tooltip {
          background: rgba(69,77,93,.9);
          border-radius: .1rem;
          color: #fff;
          display: block;
          font-size: 14px;
          max-width: 320px;
          padding: .2rem .4rem;
          position: absolute;
          text-overflow: ellipsis;
          white-space: pre;
          z-index: 300;
          visibility: hidden;
    }
    .tooltip-header {
      font-size: 1.3em;
    }
    .tooltip-key {
      font-weight: bold;
    }
    .branch path{
      stroke-width:2;
      stroke: dimgrey;
      stroke-linejoin:round;
      cursor: pointer;
      }
      .branch.hovered path{
        stroke-width:4;
        stroke: dimgrey;
      }
        .node.hovered circle{
        stroke-width:5;
        stroke: dimgrey
        }
        .node text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
          font-weight: 300;
          font-size: 0.9em;
        }
      /* .starter-template {
        padding: 40px 15px;
        text-align: left;
      } */
      .dataTables_wrapper.no-footer .dataTables_scrollBody {
        border-top: 1px solid  rgb(148, 148, 148);
        border-bottom: none;
      }
      .svg-icon {
      display: inline-flex;
      align-self: center;
      }
      h3{
          font-size: 1em;
      }
      #toTopBtn {
      position: fixed;
      bottom: 26px;
      right: 39px;
      z-index: 98;
      padding: 21px;
      background-color: ${themeColor}
      }
      .js .cd-top--fade-out {
          opacity: .5
      }
      .js .cd-top--is-visible {
          visibility: visible;
          opacity: 1
      }
      .js .cd-top {
          visibility: hidden;
          opacity: 0;
          transition: opacity .3s, visibility .3s, background-color .3s
      }
      .cd-top {
          position: fixed;
          bottom: 20px;
          bottom: var(--cd-back-to-top-margin);
          right: 20px;
          right: var(--cd-back-to-top-margin);
          display: inline-block;
          height: 40px;
          height: var(--cd-back-to-top-size);
          width: 40px;
          width: var(--cd-back-to-top-size);
          box-shadow: 0 0 10px rgba(0, 0, 0, .05) !important;
          background: url(https://res.cloudinary.com/dxfq3iotg/image/upload/v1571057658/cd-top-arrow.svg) no-repeat center 50%;
          background-color: ${themeColor};
          background-color: hsla(var(--cd-color-3-h), var(--cd-color-3-s), var(--cd-color-3-l), 0.8)
      }
      .slidecontainer {
        width: 100%;
      }
      .colourSelect {
        background: #eee;
        border-radius: 5px;
        padding: 4px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
      }
      .slider {
        -webkit-appearance: none;
        width: 100%;
        height: 15px;
        background: #d3d3d3;
        border-radius: 5px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
        -webkit-transition: .2s;
        transition: opacity .2s;
      }
      .slider:hover {
        opacity: 1; 
      }
      .slider::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 25px;
        height: 25px;
        border-radius: 50%; 
        background: ${themeColor};
        stroke: dimgrey;
        cursor: pointer;
      }
      .slider::-moz-range-thumb {
        width: 25px;
        height: 25px;
        border-radius: 50%;
        stroke: dimgrey;
        background: ${themeColor};
        cursor: pointer;
      } 
      .tree-container{
        max-height: 1000px;
        overflow: scroll;
      }
      .label{
        display: none;
      }
      .label.show{
        display: inline;
      }
      .node.hovered .label {
          display:inline;
        }
      div.sticky {
          position: -webkit-sticky; /* Safari */
          position: sticky;
          top: 0;
        }
        .searchbar {
          border-style:solid; 
          border-color: lightgrey; 
          border-radius: 5px; 
          float:right
        }
      @media print {
        .tree-container{
        max-height: none;
        overflow: visible;
        
        }
        
        .slider-block {
          display: none;
        }
        .container {
        padding-right: 1.5cm;
        padding-left: 1.5cm;
        padding-bottom: 1.5cm;
        margin: 1cm;
        min-width: 2200px;
        font-size:2.5vw;
        }
        .searchbar {
          display: none;
        }
        h3{ 
          font-size: 2.5vw;
        }
        h2 {
          font-size: 4vw;
          padding: 1cm;
        }
        h1 {
          font-size: 5vw;
        }
        .command-block {
          display: none;
        }
        pre {
          display: none;
        }
        .civet-logo {
          width: 2cm;
          height: 2cm;
        }
        .tree_svg {
          width: 1200px
        }
        .page-footer {
          display: none;
        }
        .civet-header {
          text-align: left;
        }
        .content-block, p {
        page-break-inside: avoid;
        }
      }
      @media screen and (prefers-color-scheme: dark) {
        body {
            background-color: #17141F;
            color: #F2E7DC;
            opacity: 0.95;
          }
          img {
            filter: brightness(.8) contrast(1.2);
          }
        .component {
          background-color: #2C2640;
        }
        .table-striped>tbody>tr:nth-child(odd) {
          background-color: #2C2640;
          /* border-top-color: #3B325B; */
          /* border-color: #2C2640; */
          opacity: 0.9;
        }
        .table {
          border-top: 0px;
          /* border-color: #3B325B; */
          background-color: #17141F;
      }
      .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
        border-top: none;
      }
      .accordion {
          background-color: #2C2640;
          color: #F2E7DC;
          cursor: pointer;
          padding: 12px;
          width: 100%;
          border: none;
          text-align: left;
          outline: none;
          transition: 0.4s;
        }

        .active, .accordion:hover {
          background-color: #17141F;
        }

        .accordion:after {
          content: '\002B';
          color: #F2E7DC;
          font-weight: bold;
          float: right;
          margin-left: 5px;
        }

        .active:after {
          content: "\2212";
        }

        .panel {
          padding: 0 12px;
          background-color: #2C2640;
          max-height: 0;
          overflow: hidden;
          transition: max-height 0.2s ease-out;
        }
      pre {
        background-color: #3B325B;
        color: #F2E7DC;
        border: none;
        opacity: 0.8;
      }
      .searchbar {
        background-color: #3B325B;
        color: #F2E7DC;
        border-style:none;
        opacity: 0.8;
        }
      .slider {
        background: #F2E7DC;
        stroke: #F2E7DC;
      }
      .slider::-webkit-slider-thumb {
        background: #5F9C82;
        fill: #5F9C82;
        stroke: #F2E7DC;
      }
      .slider::-moz-range-thumb {
        stroke: #F2E7DC;
        background: #5F9C82;
        fill: #5F9C82;
      } 
      .node-background{
          fill:#F2E7DC;
          stroke:#F2E7DC;
          opacity: 0.85;
      }
      .node circle{
        stroke-width:0;
        cursor:pointer;
        fill:#5F9C82;
        stroke:#F2E7DC;
        
        }
      .node circle.selected{
        stroke-width:0;
        cursor:pointer;
        fill:#E27E7E;
        stroke:#F2E7DC;
        opacity: 1;
        }
      .node rect{
        stroke-width:2;
        fill:#E27E7E;
        stroke:#F2E7DC;
      }
      .svg-tooltip {
          background: rgba(69,77,93,.9);
          color: #F2E7DC;
    }
    .branch path{
      stroke: #F2E7DC;
      opacity: 0.85;
      }
      .branch.hovered path{
        stroke:#F2E7DC;
        opacity: 1;
      }
        .node.hovered circle{
          stroke:#F2E7DC;
        opacity: 1;
        }
        .node text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
          font-weight: 300;
          font-size: 0.9em;
          color: #F2E7DC;
          fill: #F2E7DC;        
        }
        .scale-bar line {
          stroke: #F2E7DC;
        }
        .scale-bar text{
          fill: #F2E7DC;
          color:  #F2E7DC;
        }
      }
    </style>

  </head>

  <body>
    <script>
      $(document).ready(function() {
        $(window).scroll(function() {
        if ($(this).scrollTop() > 20) {
        $('#toTopBtn').fadeIn();
        } else {
        $('#toTopBtn').fadeOut();
        }
        });
        
        $('#toTopBtn').click(function() {
        $("html, body").animate({
        scrollTop: 0
        }, 400);
        return false;
        });
        });
    </script>

<script>
  function exportImageSVG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadSvg(document.querySelector(svgID), name);
      };
  };
  function exportImagePNG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadPng(document.querySelector(svgID), name);
      };
  };
</script>
    <div class="container">
      <a href="#" id="toTopBtn" class="cd-top text-replace js-cd-top cd-top--is-visible cd-top--fade-out" data-abc="true"></a>
      <div>
        <header class="civet-header">
            squirrel | 
            <small class="text-muted">Some QUIck Reconstruction to Resolve Evolutionary Links</small>
            <hr>
        </header>
        
        <h1>squirrel report
            <small class="text-muted" style="color:${themeColor}">${date}</small>
        </h1> 
        <hr>
        <h3>Best-practice phylogenetics for MPXV analysis</h3>
        <br>
        </div>
    
    <br>
    <% section_count = 1 %>
    %for clade in data_for_both_clades:
    <%data_for_report=data_for_both_clades[clade]%>
    <h2>${section_count}. Alignment for ${clade}</h2>
    <% section_count +=1 %>
      <p>
        Output alignment written to <a href="${config['epi2me_outdir']}">${config['epi2me_outdir']}</a> 
      </p>
    
    <% figure_count = 0 %>
    
    %if config["run_phylo"]:
        <h2>${section_count}. Tree file for ${clade}</h2>
        <% section_count +=1 %>
        <p>
          Output tree file written to <a href="${config['epi2me_outdir']}">${config['epi2me_outdir']}</a> 
        </p>
        %if config["run_apobec3_phylo"]:
          <h2>${section_count}. Phylo image for ${clade}</h2>
          <% section_count +=1 %>
          <br>
          <button class="accordion">Export image</button>
              <div class="panel">
              <div class="row">
                  <div class="column">
                  <button id="phylo_svg">SVG</button>
                  </div>
                  <div class="column">
                  <button id="phylo_png">PNG</button>
                  </div>
              </div>
              </div>
              <div style="width: 90%" id="phylo_div"> 
                ${data_for_report["phylo_svg_string"]}
              </div>  

          <script type="text/javascript">
          exportImageSVG("#phylo_svg","#phylo_div", "phylo_image");
          </script>
          <script type="text/javascript">
          exportImagePNG("#phylo_png","#phylo_div", "phylo_image");
          </script>
          <% figure_count +=1 %>
          <h3><strong>Figure ${figure_count}</strong> | Phylogenetic tree with APOBEC3 reconstruction for ${clade}</h3>
                  <hr>  
        %endif
    %endif
    %if config["seq_qc"]:

      <h2>${section_count}. Suggested sites to mask for ${clade}</h2>
      <% section_count +=1 %>
    
      <table class="display nowrap" id="myTable${section_count}">
        <thead>
          <tr>
           <% header = ["Name","note","present_in"] %>
            %for col in header:
              <th>${col.title().replace("_"," ")}</th>
            %endfor
          </tr>
        </thead>
        <tbody>
          <% import collections %>
          % for row in data_for_report["mask_csv"]:
            <tr>
              %for col in header:
                %if col == "present_in":
                  <td>${row[col].replace(";","<br>")}</td>
                %else:
                <td>${row[col]}</td>
                %endif
              %endfor
            </tr>
          % endfor
        </tbody>
      </table>
      <script type="text/javascript">
        $(document).ready( function () {
            var table = $('#myTable${section_count}').DataTable({
              select: {
                      style: 'multi'
                  },
              'iDisplayLength': 100,
              "paging": false,
              "border-bottom":false,
              "bInfo" : false,
              dom: 'frtip',
              buttons: ["copy","csv","print"]
            });
            table.buttons().container().appendTo( $('#tableExportID1') );
            
          } );
                  
  </script>

    %endif
  % endfor

<script>
  function exportImageSVG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadSvg(document.querySelector(svgID), name);
      };
  };
  function exportImagePNG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadPng(document.querySelector(svgID), name);
      };
  };
</script>
    <script>
        var acc = document.getElementsByClassName("accordion");
        var i;
        for (i = 0; i < acc.length; i++) {
              acc[i].addEventListener("click", function() {
                this.classList.toggle("active");
                var panel = this.nextElementSibling;
                if (panel.style.maxHeight) {
                  panel.style.maxHeight = null;
                } else {
                  panel.style.maxHeight = panel.scrollHeight*1.2 + "px";
                } 
              });
            }
    </script>

<footer class="page-footer">
  <div class="container-fluid text-right text-md-right">
    <hr>
    <div class="row">
      <div class="col-sm-1">
        <p>
        <img class="piranha-logo" src="https://raw.githubusercontent.com/artic-network/website/master/assets/images/artic-logo-small.svg" vertical-align="left" width="50" height="50"></img>
        <p>
    </div>

  <div class="col-sm-11" style="text-align: right;">
    squirrel ${config['version']} | <small class="text-muted">Some QUIck Reconstruction to Resolve Evolutionary Links</small> <br><small class="text-muted">GNU General Public License v3.0</small></div>

    <br><br>
    </p>
  </div>
</footer>
    </div>
  </body>
</html>
