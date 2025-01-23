# MiSFet
Readme.

The order of the pipeline is MIDAS2 -> Signature -> FEAST

You need to prepare several files. Initially, you need to prepare a file named “species_list.tsv” which will need to contain the MIDAS2 code of the species, which can be found in the metadata file that comes when downloading the database, which I can now see can be a problem if you’re running this a first time. But don’t fret, Bourbon has your back. I will put the metadata file from Midas2 in the repo so you can find the species codes that you need. 

Your “species_list.tsv” file. Make sure that there is one species per line and save as that name and put it in the initial folder where the Snakfile is. 

Next file you need to provide, is the sink_source.csv file that Signature uses. I will put an example of that in the repo, and please place that in the metadata folder where there will be an example of this wonderful file.

Ok, so put all your reads in the input folder or better, symlink them if you have them somewhere else so you don’t use a shit ton of space. Don’t put all the reads ever. Only reads that you want for this project obviously. Now, make sure your folder is called “input”. Technically you should be able to change this in the config, but don’t.

Now, you need, pair end reads but based on what they called, make sure to modify the config file “input_fn_pattern”. 

Speaking of the config file, there are other parameters there that you can set and change as you desire. So, change away to your heart’s desire, for that is what I am here for.

Ok, let’s start with the first step MIDAS2. Yes this uses Midas. Step by step. You start by building your initial general marker database, and you detect your species and then merge them here. Why merge species, heck if I know, but MIDAS2 website says to do it so we do it.

Ok, next, now here is where that species_list.tsv that you prepared with such hardship comes into play. We want to grab the marker that are specific to that species. Now, how MIDAS2 picks saps based on a marker genome is fun shit and we can talk about this later in detail, if you, dear user, are terribly interested because this is a terribly uninteresting discourse.

Afterwards, with these and improved and complete specific marker we go ahead to dig for all the different snps that could be present in this data.

Ok, so, next. Now, you see, not all of species are going to merge, species in your file that is because they need to be in at least two samples to merge, otherwise bye bye and the pipeline fails. Fun times. Yes, fun times.

So, now it has all your merged species and here is where MIDAS2 bids you adieu and we move on to Signature. Signature, yes. Signature pretty much takes the output from these. But wait, almost forgot, we need some name changing and unzipping and rezipping before we do this. How do we do this? Well, you don’t need to worry your pretty little heart, cause I do it for you. So, that sink_source.csv file that you used. Remember, your final file needs to be called sink_source.csv and be in the metadata folder. There is an example there currently.

So, Signature is just a python script with a bunch of variables that you can specify in the config file. Though depending on the size of your files and the variables you choose, it could take varying amounts of time.

Finally, the R script signature_to_Feast.R will give you the FEAST output.


