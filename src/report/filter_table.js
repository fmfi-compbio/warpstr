function filter_table() {
  let dropdown, table, rows, cells, country, filter;
  strand_dropdown = document.getElementById("strand_dropdown");
  table = document.getElementById("myTable");
  rows = table.getElementsByTagName("tr");
  strand_filter = strand_dropdown.value;
  for (let row of rows) {
    cells = row.getElementsByTagName("td");
    sample = cells[0] || null;
    strand = cells[1] || null;
    if ((strand_filter === "All" || !strand || (strand_filter === strand.textContent))){
      row.style.display = "";
    }
    else {
      row.style.display = "none";
    }
  }
    create_graph()
}
