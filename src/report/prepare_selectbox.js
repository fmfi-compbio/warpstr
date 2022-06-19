
select_box = document.getElementById("sample_dropdown");
var list_of_samples = [];

var table = document.getElementById("summaryTable");
var rows = table.getElementsByTagName("tr");  
for (let i = 2; i < rows.length; i++) {
  header = rows[i].getElementsByTagName("th");
  if (header.length>1){
      console.log(header[0].textContent);
      if (!list_of_samples.includes(header[0].textContent)){
            list_of_samples.push(header[0].textContent);
        }
    }
}

for (let i = 0; i < list_of_samples.length; i++) {
    select_box.options[select_box.options.length] = new Option(list_of_samples[i], list_of_samples[i]);
}
