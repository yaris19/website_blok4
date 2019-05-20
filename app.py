from flask import Flask, request, render_template
import mysql.connector

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def site():
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='fifkv@hannl-hlo-bioinformatica-mysqlsrv',
        database='fifkv',
        password='613633')
    cursor = connection.cursor()
    values = ['Header', 'Sequence', 'Accession', 'Defenition',
              'Organism_species', 'Organism_genus', 'Organism_family',
              'Ident_perc', 'Query_cov', 'E_value']
    name = request.form.get('select', False)
    limit = '10'
    if request.form.get('limit', False):
        limit = request.form.get('limit')
    search = ''
    if name == 'protein':
        search = " where defenition like '%" + request.form.get(
            'description') + "%'"
    elif name == 'organism':
        search = " where organism_species like '%" + request.form.get(
            'description') + "%' or organism_genus like '%" + request.form.get(
            'description') + "%' or organism_family like '%" + request.form.get(
            'description') + "%'"

    show_values = []
    for value in values:
        if value in request.form:
            show_values.append(value)
    select = ','.join(show_values)

    query = "select " + select + " from protein p join protein_attribute" \
                                 " pa on p.name_id = pa.name_id join " \
                                 "organism o on pa.organism_id = " \
                                 "o.organism_id join sequence s on pa.seq_id" \
                                 " = s.seq_id" + search + " limit " + limit

    print(query)
    if not show_values:
        return render_template('index.html', rows='', show_values='',
                               values=values)
    else:
        cursor.execute(query)
        rows = cursor.fetchall()
        return render_template('index.html', rows=rows, show_values=show_values,
                               values=values)


if __name__ == '__main__':
    app.run()
