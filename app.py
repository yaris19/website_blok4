from flask import Flask, request, render_template
import mysql.connector
import itertools

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('website.html')


@app.route('/database', methods=['GET', 'POST'])
def site():
    # make a connection with the database
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='fifkv@hannl-hlo-bioinformatica-mysqlsrv',
        database='fifkv',
        password='613633')
    cursor = connection.cursor()
    values_1 = ['Header', 'Definition', 'Organism_species', 'Ident_perc',
                'Accession']
    values_2 = ['Sequence', 'Organism_genus', 'Organism_family', 'Query_cov',
                'E_value', ]

    # set default limit to 10
    limit = '10'
    if request.args.get('limit', False):
        # set limit to input from user
        limit = request.args.get('limit')

    # set search to empty
    search = ''

    name = request.args.get('select', False)
    # set search to input from user
    if name == 'protein':
        search = " where definition like '%" + request.args.get(
            'description') + "%' or accession like '%" + \
                 request.args.get('description') + "%' "
    elif name == 'organism':
        search = " where organism_species like '%" + request.args.get(
            'description') + "%' or organism_genus like '%" + request.args.get(
            'description') + "%' or organism_family like '%" + request.args.get(
            'description') + "%' "

    # new list with the values that the user selected which should be shown
    # in the table
    show_values = []
    for value_1, value_2 in itertools.zip_longest(values_1[0:-1], values_2):
        if value_1 in request.args:
            # append the values that user selected to new list
            show_values.append(value_1)
        if value_2 in request.args:
            show_values.append(value_2)
    if values_1[-1] in request.args:
        show_values.append(values_1[-1])

    # seperate the values with a comma
    select = ','.join(show_values)

    # get choice from user
    choice = request.args.get('choice', False)
    # default choice is a regular select statement, so count statement is false
    count = False
    query = "select " + select + " from protein p join protein_attribute" \
                                 " pa on p.name_id = pa.name_id join " \
                                 "organism o on pa.organism_id = " \
                                 "o.organism_id join sequence s on pa.seq_id" \
                                 " = s.seq_id" + search + "limit " + limit

    if choice == 'count':
        search = " where " + select + " like '%" + request.args.get(
            'description') + "%' "
        query = "select " + select + ",count(*) from protein p join protein_" \
                                     "attribute pa on p.name_id = pa.name_id" \
                                     " join organism o on pa.organism_id = " \
                                     "o.organism_id join sequence s on " \
                                     "pa.seq_id  = s.seq_id" + search + \
                "group by " + select + " order by count(*) desc"
        # user chose count, so count statement is true
        count = True

    if not show_values:
        # query = "select header, sequence, score from sequence"
        # cursor.execute(query)
        # rows = cursor.fetchall()
        # return empty html file when site is loaded for first time
        return render_template('database.html', data='',
                               show_values=['header', 'sequence', 'score'],
                               values_1=values_1, values_2=values_2, count=None,
                               ncbi_links=None, rows='')
    else:
        # run query in database
        cursor.execute(query)
        rows = cursor.fetchall()
        ncbi_links = []
        if 'Accession' in show_values:
            for row in rows:
                link = 'https://www.ncbi.nlm.nih.gov/protein/' + row[-1]
                ncbi_links.append(link)
            data = zip(rows, ncbi_links)
            # return html file with all the variables
            return render_template('database.html', data=data,
                                   show_values=show_values,
                                   count=count, values_1=values_1,
                                   values_2=values_2, accession=True, rows=rows)
        else:
            data = rows
            return render_template('database.html', data=data,
                                   show_values=show_values,
                                   values_1=values_1, values_2=values_2,
                                   count=count, accession=False, rows=rows)


if __name__ == '__main__':
    app.run(debug=True)
