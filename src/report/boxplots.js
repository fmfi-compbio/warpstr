var trows = document.getElementById('myTable').rows;
    
for (var row = 0; row < trows.length; row++) {
    var cols = trows[row].cells;
    console.log(cols.length);
    cols[6].style.display = 'none';
    cols[7].style.display = 'none';
    cols[8].style.display = 'none';
}

function create_graph() {
    let table,rows,cells;
    var bc_data = [];
    var tr_data = [];
    var trr_data = [];
    var readnames = [];     

    strand_dropdown = document.getElementById("strand_dropdown");
    strand_filter = strand_dropdown.value;
    table = document.getElementById("myTable");
    rows = table.getElementsByTagName("tr");  
    for (let row of rows) {
      cells = row.getElementsByTagName("td");
        header = row.getElementsByTagName("th");
        if (cells.length > 3){
            strand = cells[1] || null;
            if ((strand_filter === "All" || !strand || (strand_filter === strand.textContent))){   
              bc_data.push(cells[2].textContent);
              tr_data.push(cells[3].textContent);
              trr_data.push(cells[4].textContent);
                readnames.push(header[0].textContent);
            }
        }
    }
    console.log(readnames);

    var bc_trace = {
      y: bc_data,
        text: readnames,
      type: 'box',
      name: 'TR lengths using basecalling',
      jitter: 0.7,
      pointpos: -1.8,
      marker: {
        color: 'rgb(254,74,73)'
      },
      boxpoints: 'all'
    };

    var tr_trace = {
      y: tr_data,
        text: readnames,
      type: 'box',
      name: 'TR lengths using DTW',
      jitter: 0.7,
      pointpos: -1.8,
      marker: {
        color: 'rgb(42,183,202)'
      },
      boxpoints: 'all'
    };

    var trr_trace = {
      y: trr_data,
      text: readnames,
      type: 'box',
      name: 'TR lengths using rescaled DTW',
      jitter: 0.7,
      pointpos: -1.8,
      marker: {
        color: 'rgb(254,215,102)'
      },
      boxpoints: 'all'
    };



    var data = [bc_trace, tr_trace, trr_trace];

    var layout = {
      title: 'Allele lengths'
    };


    var myPlot = document.getElementById('boxplot_div'),
        d3 = Plotly.d3,
        data = data,
        layout = {
        title: 'Allele lengths'
        };

    Plotly.newPlot('boxplot_div', data, layout);

    myPlot.on('plotly_click', function(data){
        var pts = '';
        for(var i=0; i < data.points.length; i++){
            pts = 'x = '+data.points[i].x +'\ny = '+
                data.points[i].y.toPrecision(4) + '\n\n';
            console.log(data.points[i].text);
        }
      
        alert('Closest point clicked:\n\n'+pts);
        console.log('Closest point clicked:\n\n'+pts);
    });

    }

create_graph();