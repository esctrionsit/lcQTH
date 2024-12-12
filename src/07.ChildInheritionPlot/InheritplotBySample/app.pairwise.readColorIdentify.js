const fs = require('fs');
const d3 = require("d3");
const jsdom = require("jsdom");
const { JSDOM } = jsdom;
const dom = new JSDOM();
const child_process = require('child_process');
const args = process.argv.slice(2)

// User Settings
let gIBDheader = args[0];
let gIBDline = args[1];
let Title = args[2];
let Outputfile = args[3];
let Binsize = args[4]; // Kb

// Read input
const replacer = new RegExp("\n", 'g');
gIBDheader = fs.readFileSync(gIBDheader).toString().replace(replacer, '').split("\t").slice(1);
gIBDline = fs.readFileSync(gIBDline).toString().replace(replacer, '').split("\t")[1];

// pre-process files
let Phasedata = {};
let TypeCount = [0, 0, 0, 0, 0, 0, 0];
for(let i=0; i<gIBDheader.length; i++){
	let chrom = gIBDheader[i].split(":")[0];
	if(Object.keys(Phasedata).indexOf(chrom) == -1){
		Phasedata[chrom] = [];
	}

	if(gIBDline[i] == '0'){
		Phasedata[chrom].push(0);
		TypeCount[0] += 1;
	}else if(gIBDline[i] == '1'){
		Phasedata[chrom].push(1);
		TypeCount[1] += 1;
	}else if(gIBDline[i] == '2'){
		Phasedata[chrom].push(2);
		TypeCount[2] += 1;
	}else if(gIBDline[i] == '3'){
		Phasedata[chrom].push(3);
		TypeCount[3] += 1;
	}else if(gIBDline[i] == '4'){
		Phasedata[chrom].push(4);
		TypeCount[4] += 1;
	}else if(gIBDline[i] == '5'){
		Phasedata[chrom].push(5);
		TypeCount[5] += 1;
	}else{
		Phasedata[chrom].push(6);
		TypeCount[6] += 1;
	}
}

// Load Chromosome Info
let chromlen = {};
let chrommaxlen = 0;
for(let i=0; i<Object.keys(Phasedata).length; i++){
	let chrom = Object.keys(Phasedata)[i];
	chromlen[chrom] = Phasedata[chrom].length;
	if(chrommaxlen < chromlen[chrom]){
		chrommaxlen = chromlen[chrom];
	}
}

// Load Chromosome Centromere
let centropos = fs.readFileSync("../chrcent.txt").toString().split("\n");
let chromcen = {};
let chromorder = [];
for(let i=0; i<centropos.length; i++){
	if(centropos[i].length == 0){continue;}
	let chrom = centropos[i].split("\t")[0];
	let pos = parseInt(centropos[i].split("\t")[1]);
	// chromcen[chrom] = pos
	chromcen[chrom] = pos*(1000/parseFloat(Binsize))
	chromorder.push(chrom)
}


// Plot parameters
let trackheight = 30*(1000/parseFloat(Binsize));
let binwidth = 1;
// let trackheight = 120;
// let binwidth = 0.5;
// let colors = ["#FF3D43", "#ff99ff", "#2566CB", "#dddddd", "#22bbee", "#000000"]; //, "#ccaf00"];
let colors = ["#FF3D43", "#ffffff", "#2566CB", "#dddddd", "#22bbee", "#000000"]; //, "#ccaf00"];
// let colors = ["#ffffff", "#4fb8fe", "#e01921"];
// let colors = ["#eeeeee", "#000085", "#ea3323"] // CNV color pad
let headertextwidth = 100;
let headeraxisheight = 30;
let titlefontsize = trackheight*1.5;

// init plot
let SVGcanvas = d3.select(dom.window.document.body).append("svg")
	.attr("version", "1.1")
	.attr("xmlns", "http://www.w3.org/2000/svg")
	.attr("xmlns:xlink", "http://www.w3.org/1999/xlink")
	.attr("xml:space", "preserve")
let SVGLayer = SVGcanvas.append("g");
let SVGLayerBAR = SVGLayer.append("g");
let SVGLayerCHR = SVGLayer.append("g");
let SVGLayerLegend = SVGcanvas.append("g");
let SVGheadertextLayer = SVGcanvas.append("g");
let SVGTitleLayer = SVGcanvas.append("g");
let xScale = d3.scaleLinear()
	.domain([0, chrommaxlen])
	.range([0, chrommaxlen*binwidth]);
let xAxis = d3.axisBottom(xScale)
	.tickSizeOuter(0);

// Title
SVGTitleLayer.append("text")
	.attr("x", chrommaxlen*binwidth/2)
	.attr("y", 0)	
	.attr("text-anchor", "middle")
	.attr("fill", "#000")
	.attr("font-family", "Arial")
	.style("font-size", titlefontsize+"px")
	.style("font-weight", "600")
	.text(Title);

// xAxis
let xAxisChart = SVGcanvas.append('g')
	.style("stroke-width", "2px")
	.call(xAxis);
xAxisChart.append("text")
	.attr("x", xScale(chrommaxlen) + 10)
	.attr("y", 13)	
	.attr("text-anchor", "start")
	.attr("fill", "#000")
	.attr("font-family", "Arial")
	.style("font-size", "13px")
	.text("(Mbp)");

// Legend
// SVGLayerLegend.append("rect")
// 	.attr("x", 0).attr("y", 0)
// 	.attr("width", trackheight).attr("height", trackheight)
// 	.attr("fill", colors[4]);
// SVGLayerLegend.append("text")
// 	.attr("x", trackheight + 3)
// 	.attr("y", trackheight)	
// 	.attr("text-anchor", "start")
// 	.attr("fill", "#000")
// 	.attr("font-family", "Arial")
// 	.style("font-size", trackheight + "px")
// 	.text("Shared Genomic resource Regions");
// SVGLayerLegend.append("rect")
// 	.attr("x", 300).attr("y", 0)
// 	.attr("width", trackheight).attr("height", trackheight)
// 	.attr("fill", colors[5]);
// SVGLayerLegend.append("text")
// 	.attr("x", 300 + trackheight + 3)
// 	.attr("y", trackheight)	
// 	.attr("text-anchor", "start")
// 	.attr("fill", "#000")
// 	.attr("font-family", "Arial")
// 	.style("font-size", trackheight + "px")
// 	.text("Polymorphism Hotspot Regions");

// Plot
Trackidx = 0;
for(let i=0; i<chromorder.length; i++){
	let chrom = chromorder[i];
	// Chromosome shape
	SVGLayerCHR.append("path")
		.attr("d", "M " + 0 + " " + trackheight*Trackidx + " l 0 " + trackheight + " l  " + trackheight/2 + " 0 a 1 1 0 0 1 0 -" + trackheight + " Z")
		.attr("stroke", "white").attr("stroke-width", "1px").attr("fill", "white");
	SVGLayerCHR.append("path")
		.attr("d", "M " + chromcen[chrom]*binwidth + " " + trackheight*Trackidx + " l 0 " + (trackheight) + " l  " + trackheight/2 + " 0 a 1 1 0 0 1 0 -" + trackheight + " Z")
		.attr("stroke", "white").attr("stroke-width", "1px").attr("fill", "white");
	SVGLayerCHR.append("path")
		.attr("d", "M " + chromcen[chrom]*binwidth + " " + trackheight*Trackidx + " l 0 " + (trackheight) + " l -" + (trackheight/2) + " 0 a 1 1 0 0 0 0 -" + (trackheight) + " Z")
		.attr("stroke", "white").attr("stroke-width", "1px").attr("fill", "white");
	SVGLayerCHR.append("path")
		.attr("d", "M " + chromlen[chrom]*binwidth + " " + trackheight*Trackidx + " l 0 " + (trackheight) + " l -" + (trackheight/2) + " 0 a 1 1 0 0 0 0 -" + (trackheight) + " Z")
		.attr("stroke", "white").attr("stroke-width", "1px").attr("fill", "white");
	SVGLayerCHR.append("path")
		.attr("d", "M " + trackheight/2 + " " + trackheight*Trackidx + " a 1 1 0 0 0 0 " + trackheight + " l " + (chromcen[chrom]*binwidth-trackheight) + " 0 a 1 1 0 0 0 0 -" + trackheight + " Z")
		.attr("stroke", "black").attr("stroke-width", "2px").attr("fill", "none");
	SVGLayerCHR.append("path")
		.attr("d", "M " + (chromcen[chrom]*binwidth + trackheight/2) + " " + trackheight*Trackidx + " a 1 1 0 0 0 0 " + trackheight + " l " + ((chromlen[chrom]-chromcen[chrom])*binwidth-trackheight) + " 0 a 1 1 0 0 0 0 -" + trackheight + " Z")
		.attr("stroke", "black").attr("stroke-width", "2px").attr("fill", "none");
	SVGheadertextLayer.append("text").attr("font-family", "Arial").attr("text-anchor", "end").attr("dominant-baseline", "test-after-edge").attr("font-weight", "bold").attr("font-size", trackheight+"px")
		.attr("x", headertextwidth-5).attr("y", trackheight*(Trackidx+0.8))
		.text(chrom);
	// Coverage Blocks
	{
		let sharecount = 1;
		let formertype = Phasedata[chrom][0];
		for(let j=1; j<chromlen[chrom]; j++){
			if(Phasedata[chrom][j] == formertype){
				sharecount += 1;
			}else if(sharecount > 0){
				SVGLayerBAR.append("rect")
					.attr("x", (j-sharecount)*binwidth).attr("y", trackheight*(Trackidx))
					.attr("width", binwidth*sharecount).attr("height", trackheight)
					.attr("fill", colors[formertype]);
				sharecount = 1;
				formertype = Phasedata[chrom][j];
			}
		}
		if(sharecount > 0){
			SVGLayerBAR.append("rect")
				.attr("x", (chromlen[chrom]-sharecount)*binwidth).attr("y", trackheight*(Trackidx))
				.attr("width", binwidth*sharecount).attr("height", trackheight)
				.attr("fill", colors[formertype]);
			sharecount = 0;
		}
	}
	
	Trackidx += 2;
}


SVGcanvas
	.attr("width", (chrommaxlen)*binwidth + ((chrommaxlen)*binwidth*0.2 > 40 ? (chrommaxlen)*binwidth*0.2 : 40) + headertextwidth*2)
	.attr("height", (Trackidx+1)*trackheight + ((Trackidx+1)*trackheight*0.2 > 40 ? (Trackidx+1)*trackheight*0.2 : 40) + headeraxisheight + titlefontsize + trackheight*2);
SVGLayer.attr("transform", "translate(" + (((chrommaxlen)*binwidth*0.1 > 20 ? (chrommaxlen)*binwidth*0.1 : 20) + headertextwidth) + ", " + (((Trackidx+1)*trackheight*0.1 > 20? (Trackidx+1)*trackheight*0.1 : 20) + headeraxisheight + titlefontsize + trackheight*2) + ")");
SVGheadertextLayer.attr("transform", "translate(" + (((chrommaxlen)*binwidth*0.1 > 20 ? (chrommaxlen)*binwidth*0.1 : 20)) + ", " + (((Trackidx+1)*trackheight*0.1 > 20? (Trackidx+1)*trackheight*0.1 : 20) + headeraxisheight + titlefontsize + trackheight*2) + ")");
xAxisChart.attr("transform", "translate(" + (((chrommaxlen)*binwidth*0.1 > 20 ? (chrommaxlen)*binwidth*0.1 : 20) + headertextwidth) + ", " + (((Trackidx+1)*trackheight*0.1 > 20? (Trackidx+1)*trackheight*0.1 : 20) + titlefontsize + trackheight*2) + ")");
SVGTitleLayer.attr("transform", "translate(" + (((chrommaxlen)*binwidth*0.1 > 20 ? (chrommaxlen)*binwidth*0.1 : 20) + headertextwidth) + ", " + ((Trackidx+1)*trackheight*0.1 > 20? (Trackidx+1)*trackheight*0.1 : 20) + ")");
SVGLayerLegend.attr("transform", "translate(" + (((chrommaxlen)*binwidth*0.1 > 20 ? (chrommaxlen)*binwidth*0.1 : 20) + headertextwidth) + ", " + (((Trackidx+1)*trackheight*0.1 > 20? (Trackidx+1)*trackheight*0.1 : 20) + titlefontsize) + ")");

// console.log(TypeCount)

console.error("Outputing SVG...");
fs.writeFileSync(Outputfile + '.svg', ("<?xml version=\"1.0\" encoding=\"utf-8\"?>" + dom.window.document.body.innerHTML));
console.error("Outputing SVG...Done");

console.error("Converting to PDF...");
child_process.exec('cairosvg ' + Outputfile + '.svg -o ' + Outputfile + '.pdf',function (error, stdout, stderr) {
	if (error !== null) {
		console.error('exec error: ' + error);
	}else{
		console.error("Converting to PDF...Done");
	}
});