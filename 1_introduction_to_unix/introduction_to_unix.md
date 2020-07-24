#### Introduction to Unix ####

A Unix or Unix-like operating system generally has a wide variety of build-in functions that are extremely useful for navigating between folders, exploring the contents of files, and getting various summary statistics that are useful before undertaking a bioinformatic analysis. 

Here we will explore some of the more useful Unix commands that you will find useful throughout this course and future assignments. Make sure you get comfortable using these commands, since they will make your life a lot easier later on.

If you are ever curious about how to use a command, you can type "man" right before the command name and hit enter ("man" is short for "manual" in this case, so it will give you the command's manual).  


## mkdir

The "mkdir" command will make an empty directory (folder) with a specified name. Often we may wish to organize our work into separate folders, so it's nice to be able to easily create new ones.
Let's start by creating a folder called "test_project"

mkdir test_project

When you make folders in a Unix OS it is always a good idea to keep the names simple and avoid special characters like *#%^$ etc. Special characters can sometimes cause errors when we begin to parse through folder names with python or R scripts. Also, whatever you do, DO NOT PUT SPACES IN YOUR FOLDER NAMES, as these cause lots of issues. If you need some kind of space for readability you can use an underscore instead. 
