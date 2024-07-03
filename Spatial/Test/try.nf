process myProcess {
    output:
    path '/scratch/sb14489'

    script:
    """
    python ProduceTable.py
    """
}
