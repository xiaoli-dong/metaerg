// Dimensions of sunburst.
var width = 1000;
var height = 1000;

var radius = Math.min(width, height) / 4;

// Mapping of step names to colors.
var colors = d3.scaleOrdinal(d3.schemeCategory20b);


// Total size of all segments; we set this later, after loading the data.
var totalSize = 0;

var vis = d3.select("#chart").append("svg:svg")
    .attr("width", width)
    .attr("height", height)
    .style("margin", "10px")
       .append("svg:g")
    .attr("id", "container")
    .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

var x = d3.scaleLinear().range([0, 2 * Math.PI]);
var y = d3.scaleLinear().range([0, radius]);
//var x = d3.scaleSqrt().range([0, 2 * Math.PI]);
//var y = d3.scaleSqrt().range([0, radius]);

var partition = d3.partition()

var arc = d3.arc()
    .startAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x0))); })
    .endAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x1))); })
    .innerRadius(function(d) { return Math.max(0, y(d.y0)); })
    .outerRadius(function(d) { return Math.max(0, y(d.y1)); });
// tooltip
var tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);


//d3.json("../LCM1.cds.foam.summary.json", function(error, data) {
//d3.json("tmp.json", function(error, data) {
 //   if (error) throw error;
  //  json = data;
   // alert(json.name);
   // createVisualization(json);
//});

// Main function to draw and set up the visualization, once we have the data.
function createVisualization(json) {
    
    // Bounding circle underneath the sunburst, to make it easier to detect
    // when the mouse leaves the parent g.
    var circle = vis.append("svg:circle")
	.attr("r", 1*radius)
	.style("opacity", 0)
	.attr("id", "center");
    
    var tip = vis.append("svg:text")
	.attr("id","tip")
	.style("fill","red")
	//.style("font-size",8)
	.style("text-anchor","middle")
    	.text(function(d){return json.text});


    // Turn the data into a d3 hierarchy and calculate the sums.
    var root = d3.hierarchy(json)
	.sum(function(d) { return d.size; })
	.sort(function(a, b) { return b.value - a.value; });
	//.sort(function(a, b) { return a.data.text.toLowerCase().localeCompare(b.data.text.toLowerCase()); });

    // For efficiency, filter nodes to keep only those large enough to see.
    var nodes = partition(root).descendants()
	.filter(function(d) {
            return (d.x1 - d.x0 > 0.005); // 0.005 radians = 0.29 degrees
	});
    var uniqueNames = (function(a) {
        var output = [];
        a.forEach(function(d) {
            if (output.indexOf(d.data.text) === -1) {
                output.push(d.data.text);

            }
        });
        return output;
    })(nodes);

    colors.domain(uniqueNames); // update domain colors

    vis.selectAll("path")
	.data(nodes)
    	.enter().append("g").attr("class", "node");

    var path = vis.selectAll(".node")
	.append("svg:path")
	.attr("d", arc)
	.attr("display", function(d) { return d.depth ? null : "none"; })
        .style("fill", function(d) { return colors(d.data.text); })
	.attr("fill-rule", "evenodd")
	.style("opacity", 1)
	.on("mouseover", mouseover);
 // Get total size of the tree = value of root node from partition.
    totalSize = path.datum().value;
    
    text = vis.selectAll(".node")
	.append("text")
	.attr("dy", "5px")
	.attr("dx", function(d) {
            var a = computeTextRotation(d);
	    return a >90  ? "-40px" : "40px";
	})
    	.style("text-anchor", function(d){
	    ang = (x((d.x0 + d.x1)/2) - Math.PI / 2) / Math.PI * 180;
	    if(ang > 90){
		return "end"
	    }
	    else{
		return "start"
	    }
	})
        .attr("transform", function(d) {
	    var center = arc.centroid(d);
	    return "translate(" + center[0] + "," + center[1] + ")rotate(" + computeTextRotation(d) + ")";

	})
        .text(function(d) {
	    if(!d.children){
		return d.data.text;
	    }
	    else{
		
	    }
        });
    
    // Add the mouseleave handler to the bounding circle.
    d3.select("#container").on("mouseleave", mouseleave);
   
};

// Fade all but the current sequence, and show it in the breadcrumb trail.
function mouseover(d) {

    // Fade all the segments.
    d3.selectAll("path")
	.style("opacity", 0.2);

    var sequenceArray = d.ancestors().reverse();
    sequenceArray.shift(); // remove root node from the array

    // Then highlight only those that are an ancestor of the current segment.
    vis.selectAll("path")
	.filter(function(node) {
            return (sequenceArray.indexOf(node) >= 0);
        })
	.style("opacity", 1);

    var percentage = (100 * d.value / totalSize).toPrecision(3);
    var percentageString = percentage + "%";

    if (percentage < 0.1) {
	percentageString = "< 0.1%\n";
    }

    tooltip.text(d.data.text + " " + percentageString + " count=" + d.value + " total=" + totalSize)
        .style("opacity", 0.8)
        .style("left", (d3.event.pageX) + 0 + "px")
        .style("top", (d3.event.pageY) - 0 + "px");

}

// Restore everything to full opacity when moving off the visualization.
function mouseleave(d) {

    tooltip.style("opacity", 0);

    // Transition each segment to full opacity and then reactivate it.
    d3.selectAll("path")
	.transition()
	.duration(1000)
	.style("opacity", 1)
	.on("end", function() {
            d3.select(this).on("mouseover", mouseover);
        });
}

function computeTextRotation(d) {
  ang = (x((d.x0 + d.x1)/2) - Math.PI / 2) / Math.PI * 180;
  return (ang > 90) ? 180 + ang : ang;
}

