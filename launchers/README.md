# README for launchers

This directory contains small `.desktop` launcher files which enable execution of various
programs used for HCR via Gnome desktop icons. They are simple text files, and the ones
provided here are reasonably generic.

__IMPORTANT NOTE:__ if your `HCR_configuration` repository is in a location other than
`/home/hcr/git/HCR_configuration`, you will need to edit each of the `.desktop` files and
modify the `Icon` entry to reflect your repository location.

To install the `.desktop` files, start from this `launchers` directory, then:

1. Create a `~/Desktop` directory (if needed), and copy all of the `.desktop` files from here to there.    

   ```bash
   $ mkdir -p ~/Desktop
   $ cp *.desktop ~/Desktop/
   ```
   
2. New icons should now show up on your desktop, but will be marked with a red X to indicate that they cannot be executed. For each of the new icons, right click on the icon and select "Allow launching" to make it executable.
