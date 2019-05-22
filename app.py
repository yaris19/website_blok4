from flask import Flask, request, render_template
import mysql.connector

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('website.html')


@app.route('/database.html', methods=['GET', 'POST'])
def site():
    # make a connection with the database
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='fifkv@hannl-hlo-bioinformatica-mysqlsrv',
        database='fifkv',
        password='613633')
    cursor = connection.cursor()
    values = ['Header', 'Sequence', 'Accession', 'Defenition',
              'Organism_species', 'Organism_genus', 'Organism_family',
              'Ident_perc', 'Query_cov', 'E_value']

    # set default limit to 10
    limit = '10'
    if request.form.get('limit', False):
        # set limit to input from user
        limit = request.form.get('limit')

    # set search to empty
    search = ''

    name = request.form.get('select', False)
    # set search to input from user
    if name == 'protein':
        search = " where defenition like '%" + request.form.get(
            'description') + "%'"
    elif name == 'organism':
        search = " where organism_species like '%" + request.form.get(
            'description') + "%' or organism_genus like '%" + request.form.get(
            'description') + "%' or organism_family like '%" + request.form.get(
            'description') + "%'"

    # new list with the values that should be shown in the table
    show_values = []
    for value in values:
        if value in request.form:
            # append the values that user selected to new list
            show_values.append(value)
    # seperate the values with a comma
    select = ','.join(show_values)

    # get choice from user
    choice = request.form.get('choice', False)
    # default choice is a regular select statement, so count statement is false
    count = False
    query = "select " + select + " from protein p join protein_attribute" \
                                 " pa on p.name_id = pa.name_id join " \
                                 "organism o on pa.organism_id = " \
                                 "o.organism_id join sequence s on pa.seq_id" \
                                 " = s.seq_id" + search + " limit " + limit

    if choice == 'count':
        query = "select " + select + ",count(*) from protein p join protein_" \
                                     "attribute pa on p.name_id = pa.name_id" \
                                     " join organism o on pa.organism_id = " \
                                     "o.organism_id join sequence s on " \
                                     "pa.seq_id  = s.seq_id group by " + \
                select + " order by count(*) desc"
        # user chose count, so count statement is true
        count = True

    if not show_values:
        # return empty html file when site is loaded for first time
        return render_template('database.html', rows='', show_values='',
                               values=values, count='')
    else:
        # run query in database
        cursor.execute(query)
        rows = cursor.fetchall()
        # return html file with all the variables
        return render_template('database.html', rows=rows,
                               show_values=show_values,
                               values=values, count=count)


if __name__ == '__main__':
    app.run(debug=True)
