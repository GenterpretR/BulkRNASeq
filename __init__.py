from flask import Flask
import compute_result

app = Flask(__name__)

@app.route('/BulkRNASeq/compute_result')
def result():
    return compute_result.compute_result()
# @app.route('/gexp_to_testData')  # ask where to put the .csv file and smiles list
# def gexp_smile_to_testSet(gene_expression, smiles_list):  #GSVA Gene expression dataframe and SMILES in form of a list

if __name__ == "__main__":
    app.run(host='127.0.0.1', port=8080, debug=True)
