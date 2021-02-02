# ACM Research Coding Challenge (Spring 2021)

# Addons used:
    Python Imaging Library (PIL)
    Biopython
    ReportLab

# Process
The way I tackle a problem is, on first glance, I break it up into problems that need to be solved. 
First I skimmed through and thought, "What is a gen bank file and how do I use and manipulate it? Are there any languages that are more suited to this than other languages? What is a circular genome map and how can I pull that information out of the gen bank file?"
After asking these questions to myself, I decided to do a bit of research on what these terms actually were because I was completely in the dark about anything biology related. Knowing I could use google to understand whatever I needed to, I went right to understand what exactly a Genome Map looked like and what the components of the map look like. Next i went through and tried to understand what a GenBank file is and how to manipulate it and while going through and researching I found this one really useful article,  https://bit.ly/3jaaSJV, which led me to a useful extension in python, Biopython.
When I learned about biopython, I thought to myself, "well... I dont know how to use python but this does it for me!" So I spent two weeks learning how to use python to the point where I could use BioPython without having to struggle. I also began reading the BioPython documentation, https://biopython.org/wiki/Documentation, and found that it had all the necessary tools I would need to understand and output the genbank files.
I also understood that I would need to be able to manipulate the png file itself, so I knew I would need to find some stuff that helped me manipulate the file itself. So i looked into several things and found that Python Image Library would allow me to open the file and use it to how I would need it, and ReportLab to give me the ability to color the labels that I would want added onto the file itself. 
So once I had everything ready, I opened pycharm and used my terminal to install each of the necessary addons, I started understand what I needed to write out, and how I needed to write them out. Using the tools from biopython I was able to pull everything necessary from the genbank file and convert it using two for loops, one to pull out the genome information to make the file, and another to add labels to specific gene sequences. The BioPython algorithms I used allow me to open the file and look for specific labels that come with every genbank file, and let me categorize the information pulled however I so please. I did this by taking a for loop and going through each feature within the file and checking if the title was "gene", if it wasn't the loop moved on, if it was it assigned a color to that gene set and kept going. The 2nd for loop took two Gene Sequences and their respective names and searched for them throughout the Gene Sequence of the virus and applied a label to the respective area they fell in. 
The last part of my program adds a title to the file, which is where I used PIL's features to open and write onto a png file. I did so because I wanted to be able to add some depth to my file, just to make it feel like a completed project rather than something I was required to do.
Realistically, I should have done this in a language that I knew beforehand, like Java or C++, but I knew that Python would be a good language for me to pick up and to learn in general, so I am glad I spent the time and effort into learning something new.
