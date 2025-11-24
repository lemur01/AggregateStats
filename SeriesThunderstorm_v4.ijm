#@File (label="Input file", style="open") filename
#@String (label="Mode", choices={"Same","8-bit","16-bit","32-bit"}, description="Select the bit-depth of the image") mode
#@File (label="Output folder", style="directory") imageFolder
#@String (label="Format",choices={"TIFF","PNG","JPEG"},style="list") format
#@String(label="Tag", value="", description="Add a tag to the filename before the extension when saving the image") tag
#@boolean (label="Dummy run",value=true) dummy

/* Convert all series in a file to TIFs in a folder
 *  
 *  You can batch process files using the Batch function in the script editor
 *  
 *  Using the dummy mode enable to inspect image size as well and save the result in a cvs file
 *  
 * Jerome Boulanger for Marta 2021, updated for Nathan
 * Modified by Leila for Shekhar
 */
 
run("Close All");
setBatchMode("hide");
name = File.getNameWithoutExtension(filename);
ext = getNewFileExtension(format);

run("Bio-Formats Macro Extensions");
print("Input file:" + filename);
Ext.setId(filename);
Ext.getSeriesCount(seriesCount);
print("File contains " + seriesCount + " series");

nR0 = nResults;
for (s = 5; s <= 7; s++) if (s!=6) {
	
	
	str="open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+"";	
	
	Ext.setSeries(s-1);
	Ext.getSeriesName(seriesName);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);
	Ext.getSizeZ(sizeZ);
	Ext.getSizeC(sizeC);
	Ext.getSizeT(sizeT);
		
	print("Loading series " + s + "/" + seriesCount);			
	oname = imageFolder + File.separator + name + "_series_" + IJ.pad(s,4) + tag + ext;
	onamecsv = imageFolder + File.separator + name + "_series_" + IJ.pad(s,4) + tag + ".csv";

		
	if (!dummy) {
		run("Bio-Formats Importer", str); 
			
		idc = getImageID();								
		print("Saving series "+ s +" to "+ oname);	
		saveAs(format,oname);	
		selectImage(idc);
		id1 = processImage(onamecsv);

		
	}

	print(onamecsv);

	nR0 = nR0 + 1;
	run("Close All");
}
Ext.close();
setBatchMode("exit and display");


function getNewFileExtension(format) {
	f = newArray("TIFF","PNG","JPEG");
	e = newArray(".tif",".png",".jpg");
	for (k = 0; k < f.length; k++) {
		if(matches(format, f[k])){
			k0 = k;
			return e[k];
		}
	}
	return ".tif";
}

function processImage(onamecsv) {
	run("Camera setup", "offset=414.0 isemgain=false photons2adu=3.6 pixelsize=20.0");
	id0 = getImageID();
	run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=1.6 fitradius=3 method=[Weighted Least squares] full_image_fitting=false mfaenabled=true keep_same_intensity=false nmax=5 fixed_intensity=true expected_intensity=500:2500 pvalue=1.0E-6 renderer=[Averaged shifted histograms] magnification=5.0 colorizez=false threed=false shifts=2 repaint=50");
	run("Export results", "filepath=["+onamecsv+"] fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=true offset=true saveprotocol=true x=true y=true bkgstd=true id=false uncertainty=true frame=false");
	
	return getImageID();
}
