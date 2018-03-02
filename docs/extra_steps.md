# Introduction
After installation of R, additional steps might be require for MAC and Window users. Below are the instructions

# MAC Users
1. Download and install the latest [XQuartz](https://www.xquartz.org)
   - This is because MAC no longer ship the X11 package which is required by R to perform plotting.
2. Run `xcode-select --install` on your terminal
   -  This will install the required zlib package on your system, which is required by PRSice (for decompressing bgen files)

# Window Users
As installation of R does not automatically add it to the system path, one will need to type the full path of the R.exe and Rscript.exe in order to use PRSice. To avoid this complication, we can manually add the folder containing the R binary to the system path:

## For Windows 8 and 10
1. In Search, search for and then select: **_System (Control Panel)_**
2. Click the **_Advanced system settings_** link
3. Click the **_Advanced_** tab
4. Click _**Environment Variables**_
5. Under _**System Variables**_, select **_path_** (If you cannot find path, you can click new to make it)
6. Click _**Edit**_
7. Click _**Browse**_ and select the location of the executable of R.
If you use the default installation path, you can add `C:\Program Files\R\R-3.3.2\bin`, where (eg.) _**3.3.2**_ is the version number. Some installation might also have a i384 and x64 version and either one of those will work.

## For Windows 7
1. From the desktop, right click the Computer icon.
2. Choose _**Properties**_ from the context menu.
3. Click the _**Advanced system settings**_ link.
4. Click _**Environment Variables**_. In the section _**System Variables**_, find the _**PATH**_ environment variable and select it. Click Edit.
5. In the _**Edit System Variable**_ (or _**New System Variable**_) window, add the location of the executable of R.
If you use the default installation path, you can add `C:\Program Files\R\R-3.3.2\bin`, where (eg.) _**3.3.2**_ is the version number. Some installation might also have a i384 and x64 version and either one of those will work.
