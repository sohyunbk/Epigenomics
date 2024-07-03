process myProcess {
    output:
    path '/scratch/sb14489/Test.txt'

    script:
    """
    echo "Running ProduceTable.py"
    python ProduceTable.py
    echo "Script execution finished"
    ls -l /scratch/sb14489/
    """
}
