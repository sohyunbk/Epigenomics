process myProcess {
    output:
    path '/scratch/sb14489/Test.txt'

    script:
    """
    python ProduceTable.py
    """
}
