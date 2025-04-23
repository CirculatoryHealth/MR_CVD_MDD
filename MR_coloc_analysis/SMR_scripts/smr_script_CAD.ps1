
# This script runs SMR analysis for 18 pre-identified probes with corresponding p-value thresholds for instrument selection

# Files to have/Paths to adjust:
# - ld reference files
# - gwas file (here DEP) in COJO format
# - probe file list named probe + counter + .flist (here contain only one probe each)
# - smr software for your OS (no such thing as 'obvious')
#
# Other prep required:
# - decide on a pvalue threshold for each probe file list (here calculated by eight_regions.R)
# - create .flist files and the .esd files they refer to


# Define the range of files to process
$fileNumbers = 1..18  # Adjust this range as needed
# Define the p-values (identified in Rscript eight_regions.R)
$pvalues = @(1.05e-06, 8.30e-04, 6.25e-04, 1.14e-03, 8.04e-05, 6.65e-05, 2.42e-05, 1.23e-04, 4.47e-04, 1.96e-03, 7.54e-04, 3.55e-05, 3.55e-03, 2.49e-04, 8.88e-03, 2.54e-03, 8.15e-04, 7.14e-03)

# Loop through each file
foreach ($i in $fileNumbers) {
    # Define file names
    # $inputfile = "probe$i.flist"
    $outputfile = "../eQTL_DEP/probe$i"
    $smrfile = "probe$i.smr"
    $logfile = "probe$i.log"  # Log file for the current probe

    # Select p-value threshold
    $pvalueThreshold = $pvalues[$i - 1]
    Write-Output "Processing probe$i with p-value threshold: $pvalueThreshold"
    Add-Content $logfile "Processing probe$i with p-value threshold: $pvalueThreshold"

    # Make .besd file
    #Add-Content $logfile "Starting .besd file creation for $inputfile at $(Get-Date)"
    #./smr-1.3.1-win --eqtl-flist `"$inputfile`" --make-besd --out `"$outputfile`"
    #if ($LASTEXITCODE -ne 0) {
    #    $errorMessage = "Error creating .besd file for $inputfile at $(Get-Date)"
    #   Write-Host $errorMessage
    #    Add-Content $logfile $errorMessage
    #    continue  # Skip to the next iteration if there's an error
    #} else {
    #    Add-Content $logfile "Successfully created .besd file for $inputfile at $(Get-Date)"
    #}

    # Run SMR
    Add-Content $logfile "Starting SMR for $smrfile at $(Get-Date)"
    ./smr-1.3.1-win --bfile ../../../ld-files/g1000_eur --gwas-summary ../../Data/CAD/CAD_COJO.txt --beqtl-summary `"$outputfile`" --peqtl-smr `"$pvalueThreshold`" --out `"$smrfile`"
    if ($LASTEXITCODE -ne 0) {
        $errorMessage = "Error running SMR for $smrfile at $(Get-Date)"
        Write-Host $errorMessage
        Add-Content $logfile $errorMessage
        continue  # Skip to the next iteration if there's an error
    } else {
        Add-Content $logfile "Successfully completed SMR for $smrfile at $(Get-Date)"
    }
}

